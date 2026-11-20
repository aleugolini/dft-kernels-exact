/*
 * fs_search.c -- exact free-scalar linear-circuit size of small DFT / WHT
 * kernels (or any matrix) over GF(p).
 *
 * Model: gates z = a*u + b*v, a,b in GF(p) free, fan-in 2.  A circuit is a set W
 * of auxiliary points of PG(r-1,p) such that inputs + W + targets can be ordered
 * with each new point on the line through two earlier ones; L = r + min|W|
 * (targets distinct, not unit vectors, no repeated rows).
 *
 * Search: iterative deepening on |W|; targets added greedily when constructible
 * (monotone, hence safe); W enumerated in canonical "smallest constructible id
 * first" order; symmetry of the DFT (affine index maps + modulations) used on the
 * first auxiliary point; on the last level only points that make some missing
 * target constructible are tried.  Optional work split for parallel runs.
 *
 * Fields: GF(q) with q = p (prime) or q = p^2 (elements a + b*sqrt(n), n the smallest
 * non-residue mod p, encoded as a + p*b).
 *
 * build: cc -O3 -std=gnu99 -o fs_search fs_search.c   (add -march=native when built where it runs)
 * usage: fs_search r q wmin wmax [dft|wht|file:PATH] [split k K] [splitlev L] [nosym] [maxnodes N]
 *        file format: m, then m rows of r field codes (targets; inputs are unit vectors)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define MAXR 8
#define MAXD 24
#define MAXT 24
#define MAXQ 128

static int r, p, NP, NW;          /* dimension, field size q, #points, #words  */
static int ADD[MAXQ][MAXQ], MUL[MAXQ][MAXQ];
static long long maxnodes = 0; static int aborted = 0;
static uint8_t *pts;              /* NP x r normalized vectors                 */
static int64_t offs[MAXR + 1];    /* index offsets by leading position        */
static int inv[MAXQ];
static uint64_t *Tmask;           /* targets                                  */
static int targets[MAXT], ntarg;
static int inputs[MAXR];
static long long nodes;
static int split_k = 0, split_K = 1, split_lev = 1; static long long split_ctr = 0;
static int use_sym = 1;
static int print_all = 0, nprinted = 0;
static int G[2048][MAXR], Gd[2048][MAXR], nG = 0; /* sigma(n), d_n */

static int idx_of(const int *v) {           /* v reduced mod p, not all zero */
    int i = 0; while (i < r && v[i] == 0) i++;
    int s = inv[v[i]]; int64_t code = 0, mul = 1;
    for (int j = i + 1; j < r; j++) { code += (int64_t)MUL[v[j]][s] * mul; mul *= p; }
    return (int)(offs[i] + code);
}

static void build_points(void) {
    NP = 0; offs[0] = 0;
    for (int i = 0; i < r; i++) { int64_t c = 1; for (int j = i + 1; j < r; j++) c *= p; offs[i + 1] = offs[i] + c; }
    NP = (int)offs[r]; NW = (NP + 63) / 64;
    pts = malloc((size_t)NP * r);
    for (int i = 0; i < r; i++) {
        int64_t cnt = offs[i + 1] - offs[i];
        for (int64_t c = 0; c < cnt; c++) {
            uint8_t *q = pts + (offs[i] + c) * r; int64_t t = c;
            for (int j = 0; j < r; j++) q[j] = 0;
            q[i] = 1;
            for (int j = i + 1; j < r; j++) { q[j] = t % p; t /= p; }
        }
    }
}

typedef struct { int P[64]; int np; uint64_t *pm, *reach; } State;

static void st_alloc(State *s) { s->pm = calloc(NW, 8); s->reach = calloc(NW, 8); s->np = 0; }
static void st_copy(State *d, const State *s) {
    memcpy(d->P, s->P, sizeof(int) * s->np); d->np = s->np;
    memcpy(d->pm, s->pm, 8 * NW); memcpy(d->reach, s->reach, 8 * NW);
}
#define BIT(m, i) (((m)[(i) >> 6] >> ((i) & 63)) & 1)
#define SET(m, i) ((m)[(i) >> 6] |= 1ULL << ((i) & 63))

static void line_or(uint64_t *m, int u, int c) {     /* points of line(u,c) */
    const uint8_t *U = pts + (size_t)u * r, *Cv = pts + (size_t)c * r;
    int v[MAXR];
    SET(m, u);
    for (int a = 0; a < p; a++) {
        for (int j = 0; j < r; j++) v[j] = ADD[MUL[a][U[j]]][Cv[j]];
        SET(m, idx_of(v));
    }
}

static void add_point(State *s, int c) {
    for (int i = 0; i < s->np; i++) line_or(s->reach, s->P[i], c);
    s->P[s->np++] = c; SET(s->pm, c);
}

static void closure(State *s) {
    for (;;) {
        int added = 0;
        for (int t = 0; t < ntarg; t++) {
            int x = targets[t];
            if (!BIT(s->pm, x) && BIT(s->reach, x)) { add_point(s, x); added = 1; }
        }
        if (!added) return;
    }
}

static int all_targets(const State *s) {
    for (int t = 0; t < ntarg; t++) if (!BIT(s->pm, targets[t])) return 0;
    return 1;
}

static int sym_image(int g, int c) {
    const uint8_t *v = pts + (size_t)c * r; int w[MAXR];
    for (int n = 0; n < r; n++) w[G[g][n]] = MUL[Gd[g][n]][v[n]];
    return idx_of(w);
}
static int is_rep(int c) {
    for (int g = 0; g < nG; g++) if (sym_image(g, c) < c) return 0;
    return 1;
}

static State stk[MAXD];
static uint64_t *tmpmask;
static int sol[MAXD];

static int dfs(int lev, int wleft, int last, const uint64_t *Rprev) {
    State *s = &stk[lev];
    nodes++;
    if (maxnodes && nodes > maxnodes) { aborted = 1; return 0; }
    if (all_targets(s)) {
        if (!print_all) return 1;
        printf("  W:");
        for (int l = 0; l < lev; l++) {
            const uint8_t *q = pts + (size_t)sol[l] * r; printf(" (");
            for (int j = 0; j < r; j++) printf("%d%s", q[j], j < r - 1 ? "," : ")");
        }
        printf("\n");
        return ++nprinted >= print_all;
    }
    if (wleft == 0) return 0;
    uint64_t *filt = NULL;
    if (wleft == 1) {                          /* last W point must unlock a target */
        filt = tmpmask; memset(filt, 0, 8 * NW);
        for (int t = 0; t < ntarg; t++) if (!BIT(s->pm, targets[t]))
            for (int i = 0; i < s->np; i++) line_or(filt, targets[t], s->P[i]);
    }
    for (int w = 0; w < NW; w++) {
        uint64_t cand = s->reach[w] & ~s->pm[w] & ~Tmask[w];
        if (filt) cand &= filt[w];
        while (cand) {
            int b = __builtin_ctzll(cand); cand &= cand - 1;
            int c = w * 64 + b;
            if (c >= NP) break;
            if (Rprev && c < last && BIT(Rprev, c)) continue;
            if (lev == 0 && use_sym && !is_rep(c)) continue;
            if (lev == split_lev) { if ((split_ctr++ % split_K) != split_k) continue; }
            State *ch = &stk[lev + 1];
            st_copy(ch, s); add_point(ch, c); closure(ch);
            sol[lev] = c;
            /* after a symmetry-reduced first choice, ordering restarts */
            int ok = (lev == 0 && use_sym) ? dfs(lev + 1, wleft - 1, -1, NULL)
                                           : dfs(lev + 1, wleft - 1, c, s->reach);
            if (ok) return 1;
            if (aborted) return 0;
        }
    }
    return 0;
}

static int powmod(int a, int e) { int x = 1; while (e) { if (e & 1) x = MUL[x][a]; a = MUL[a][a]; e >>= 1; } return x; }

static void build_field(int q) {
    int pp = 0;
    for (int c = 2; c <= q; c++) if (q % c == 0) { pp = c; break; }
    int k = 0, t = q; while (t % pp == 0) { t /= pp; k++; }
    if (t != 1 || k > 2 || q > MAXQ) { fprintf(stderr, "q must be p or p^2 <= %d\n", MAXQ); exit(1); }
    if (k == 1) {
        for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) { ADD[a][b] = (a + b) % q; MUL[a][b] = a * b % q; }
    } else {
        int n = 0;                                   /* smallest non-residue mod pp */
        if (pp == 2) { fprintf(stderr, "char 2 not supported\n"); exit(1); }
        for (int c = 2; c < pp && !n; c++) { int sq = 0; for (int y = 1; y < pp; y++) if (y * y % pp == c) sq = 1; if (!sq) n = c; }
        for (int x = 0; x < q; x++) for (int y = 0; y < q; y++) {
            int a1 = x % pp, b1 = x / pp, a2 = y % pp, b2 = y / pp;
            ADD[x][y] = (a1 + a2) % pp + pp * ((b1 + b2) % pp);
            int re = (a1 * a2 + n * b1 * b2) % pp, im = (a1 * b2 + a2 * b1) % pp;
            MUL[x][y] = re + pp * im;
        }
    }
    for (int a = 1; a < q; a++) for (int b = 1; b < q; b++) if (MUL[a][b] == 1) inv[a] = b;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s r q wmin wmax [dft|wht|file:PATH] [split k K] [splitlev L] [nosym] [maxnodes N]\n", argv[0]); return 1; }
    r = atoi(argv[1]); p = atoi(argv[2]); int wmin = atoi(argv[3]), wmax = atoi(argv[4]);
    const char *kind = argc > 5 ? argv[5] : "dft";
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "split") && i + 2 < argc) { split_k = atoi(argv[i + 1]); split_K = atoi(argv[i + 2]); i += 2; }
        else if (!strcmp(argv[i], "splitlev") && i + 1 < argc) { split_lev = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "nosym")) use_sym = 0;
        else if (!strcmp(argv[i], "maxnodes") && i + 1 < argc) { maxnodes = atoll(argv[i + 1]); i++; }
        else if (!strcmp(argv[i], "all") && i + 1 < argc) { print_all = atoi(argv[i + 1]); i++; }
    }
    build_field(p);
    build_points();
    Tmask = calloc(NW, 8); tmpmask = calloc(NW, 8);
    int rows[MAXT][MAXR], omega = 0;
    ntarg = r;
    if (!strcmp(kind, "dft")) {
        for (int g = 2; g < p && !omega; g++) {
            if (powmod(g, r) != 1) continue;
            int prim = 1; for (int k = 1; k < r; k++) if (powmod(g, k) == 1) prim = 0;
            if (prim) omega = g;
        }
        if (!omega) { fprintf(stderr, "no element of order %d in GF(%d)\n", r, p); return 1; }
        for (int k = 0; k < r; k++) for (int n = 0; n < r; n++) rows[k][n] = powmod(omega, k * n);
    } else if (!strcmp(kind, "wht")) {
        int m1 = 0; for (int x = 0; x < p; x++) if (ADD[x][1] == 0) m1 = x;
        for (int k = 0; k < r; k++) for (int n = 0; n < r; n++) rows[k][n] = (__builtin_popcount(k & n) & 1) ? m1 : 1;
        use_sym = 0;
    } else if (!strncmp(kind, "file:", 5)) {
        FILE *f = fopen(kind + 5, "r"); if (!f) { perror("file"); return 1; }
        if (fscanf(f, "%d", &ntarg) != 1 || ntarg > MAXT) return 1;
        for (int k = 0; k < ntarg; k++) for (int n = 0; n < r; n++) { if (fscanf(f, "%d", &rows[k][n]) != 1) return 1; }
        fclose(f); use_sym = 0;
    } else { fprintf(stderr, "unknown kind\n"); return 1; }
    { int nt = 0;                      /* dedupe; drop targets that are inputs */
      for (int k = 0; k < ntarg; k++) {
          int id = idx_of(rows[k]), dup = 0, nz = 0;
          for (int n = 0; n < r; n++) nz += rows[k][n] != 0;
          for (int t = 0; t < nt; t++) if (targets[t] == id) dup = 1;
          if (!dup && nz > 1) { targets[nt++] = id; SET(Tmask, id); }
      }
      ntarg = nt; }
    for (int i = 0; i < r; i++) { int v[MAXR] = {0}; v[i] = 1; inputs[i] = idx_of(v); }

    if (use_sym) {               /* x -> sigma(n)=a n + b, d_n = omega^(c n) */
        for (int a = 1; a < r; a++) {
            int g = 1; for (int t = 2; t <= a; t++) if (a % t == 0 && r % t == 0) g = 0;
            if (!g) continue;
            for (int b = 0; b < r; b++) for (int c = 0; c < r; c++) {
                for (int n = 0; n < r; n++) { G[nG][n] = (a * n + b) % r; Gd[nG][n] = powmod(omega, c * n); }
                nG++;
            }
        }
        /* sanity: every group element maps the target set onto itself */
        for (int g = 0; g < nG; g++) for (int k = 0; k < ntarg; k++) {
            int im = sym_image(g, targets[k]), ok = 0;
            for (int t = 0; t < ntarg; t++) if (targets[t] == im) ok = 1;
            if (!ok) { fprintf(stderr, "symmetry check failed\n"); return 1; }
        }
    }
    for (int i = 0; i < MAXD; i++) st_alloc(&stk[i]);
    for (int i = 0; i < r; i++) add_point(&stk[0], inputs[i]);
    closure(&stk[0]);
    printf("r=%d p=%d kind=%s omega=%d points=%d |G|=%d split %d/%d level %d\n", r, p, kind, omega, NP, nG, split_k, split_K, split_lev);
    fflush(stdout);
    for (int w = wmin; w <= wmax; w++) {
        clock_t t0 = clock(); nodes = 0; split_ctr = 0;
        int ok = dfs(0, w, -1, NULL);
        printf("  |W|<=%d (L<=%d): %s  nodes=%lld  %.1fs\n", w, ntarg + w, ok ? "FOUND" : (aborted ? "ABORT" : "none"), nodes,
               (double)(clock() - t0) / CLOCKS_PER_SEC);
        if (ok && !print_all) {
            /* the successful branch is the chain of sol[] down to the leaf */
            printf("  W:");
            for (int l = 0; l < MAXD; l++) {
                if (l >= w) break;
                const uint8_t *q = pts + (size_t)sol[l] * r; printf(" (");
                for (int j = 0; j < r; j++) printf("%d%s", q[j], j < r - 1 ? "," : ")");
            }
            printf("\n"); fflush(stdout); return 0;
        }
        fflush(stdout);
    }
    return 0;
}
