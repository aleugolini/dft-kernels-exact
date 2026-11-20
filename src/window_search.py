#!/usr/bin/env python3
"""Exact local improvement of free-scalar circuits for the N-point DFT over GF(q).

1. Build several N-gate*log2N circuits for F_N (Cooley-Tukey compositions of
   2-point and 4-point kernels; the 4-point kernels include non-FFT 8-gate
   circuits found by fs_search).  Twiddles are absorbed into gate coefficients.
2. Enumerate convex windows S (connected gate sets, no path leaving and
   re-entering S) with at most DMAX independent window inputs.
3. For each window, express its outputs (gates used outside S or final outputs)
   in the coordinates of its inputs and ask fs_search whether they can be
   computed with |S|-1 gates.  A hit is a circuit for F_N with one gate fewer.

usage: python3 window_search.py q [N=8] [DMAX=5] [SMAX=10] [MAXNODES=2e7] [WORKERS=12]
N=5 starts instead from Lev's 15-gate and Winograd's 17-gate circuits (needs 5 | q-1).
DUMP=dir: write the window problems (fs_search input files + manifest.txt) and exit.
Field encoding identical to fs_search.c (q = p or p^2).
"""
import itertools
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- GF(q)
class GF:
    def __init__(self, q):
        self.q = q
        p = next(c for c in range(2, q + 1) if q % c == 0)
        k, t = 0, q
        while t % p == 0:
            t //= p; k += 1
        assert t == 1 and k <= 2
        self.p = p
        if k == 1:
            self.add = [[(a + b) % q for b in range(q)] for a in range(q)]
            self.mul = [[(a * b) % q for b in range(q)] for a in range(q)]
        else:
            n = next(c for c in range(2, p) if all(y * y % p != c for y in range(1, p)))
            self.add = [[0] * q for _ in range(q)]
            self.mul = [[0] * q for _ in range(q)]
            for x in range(q):
                for y in range(q):
                    a1, b1, a2, b2 = x % p, x // p, y % p, y // p
                    self.add[x][y] = (a1 + a2) % p + p * ((b1 + b2) % p)
                    self.mul[x][y] = (a1 * a2 + n * b1 * b2) % p + p * ((a1 * b2 + a2 * b1) % p)
        self.neg = [next(y for y in range(q) if self.add[x][y] == 0) for x in range(q)]
        self.inv = [0] + [next(y for y in range(q) if self.mul[x][y] == 1) for x in range(1, q)]

    def pw(self, a, e):
        x = 1
        for _ in range(e):
            x = self.mul[x][a]
        return x

    def root(self, N):
        for g in range(2, self.q):
            if self.pw(g, N) == 1 and all(self.pw(g, k) != 1 for k in range(1, N)):
                return g
        raise ValueError("no primitive root of order %d" % N)

    def lin(self, a, u, b, v):
        return tuple(self.add[self.mul[a][x]][self.mul[b][y]] for x, y in zip(u, v))

    def normalize(self, v):
        for c in v:
            if c:
                s = self.inv[c]
                return tuple(self.mul[s][x] for x in v)
        return None

    def solve(self, B, t):
        """coordinates c with sum_i c_i B[i] = t (B rows independent), else None."""
        m, n = len(B), len(t)
        # augmented system: columns = B rows, rhs = t
        A = [[B[i][j] for i in range(m)] + [t[j]] for j in range(n)]
        row, piv = 0, []
        for col in range(m):
            pr = next((i for i in range(row, n) if A[i][col]), None)
            if pr is None:
                return None
            A[row], A[pr] = A[pr], A[row]
            s = self.inv[A[row][col]]
            A[row] = [self.mul[s][x] for x in A[row]]
            for i in range(n):
                if i != row and A[i][col]:
                    f = self.neg[A[i][col]]
                    A[i] = [self.add[x][self.mul[f][y]] for x, y in zip(A[i], A[row])]
            piv.append(col); row += 1
        if any(A[i][m] for i in range(row, n)):
            return None
        return [A[i][m] for i in range(m)]

    def rank(self, rows):
        A = [list(r) for r in rows]
        rk, ncol = 0, len(A[0]) if A else 0
        for col in range(ncol):
            pr = next((i for i in range(rk, len(A)) if A[i][col]), None)
            if pr is None:
                continue
            A[rk], A[pr] = A[pr], A[rk]
            s = self.inv[A[rk][col]]
            A[rk] = [self.mul[s][x] for x in A[rk]]
            for i in range(len(A)):
                if i != rk and A[i][col]:
                    f = self.neg[A[i][col]]
                    A[i] = [self.add[x][self.mul[f][y]] for x, y in zip(A[i], A[rk])]
            rk += 1
        return rk


# ---------------------------------------------------------------- circuits
# A circuit: n inputs, gates [(a, u, b, v)], outputs [node ids] in DFT order.
def kernel2(F):
    return (2, [(1, 0, 1, 1), (1, 0, F.neg[1], 1)], [2, 3])


def compose(F, N1, C1, N2, C2, w):
    """DIT Cooley-Tukey: F_N, N = N1 N2, from circuits C1 (length N1) and C2 (N2).
    X[k1 + N1 k2] = sum_n2 W_N2^{n2 k2} W_N^{n2 k1} U_k1[n2],
    U_k1[n2] = sum_n1 x[n2 + N2 n1] W_N1^{n1 k1};  w = primitive N-th root."""
    N = N1 * N2
    gates, U = [], {}
    nxt = N
    for n2 in range(N2):                      # inner N1-point DFTs
        ids = {i: n2 + N2 * i for i in range(N1)}
        for g, (a, u, b, v) in enumerate(C1[1]):
            gates.append((a, ids[u], b, ids[v])); ids[N1 + g] = nxt; nxt += 1
        for k1, o in enumerate(C1[2]):
            U[(k1, n2)] = ids[o]
    out = [None] * N
    for k1 in range(N1):                      # outer N2-point DFTs, twiddles absorbed
        ids, sc = {}, {}
        for n2 in range(N2):
            ids[n2] = U[(k1, n2)]; sc[n2] = F.pw(w, n2 * k1)
        for g, (a, u, b, v) in enumerate(C2[1]):
            aa = F.mul[a][sc.get(u, 1)]; bb = F.mul[b][sc.get(v, 1)]
            gates.append((aa, ids[u], bb, ids[v])); ids[N2 + g] = nxt; sc[N2 + g] = 1; nxt += 1
        for k2, o in enumerate(C2[2]):
            out[k1 + N1 * k2] = ids[o]
    return (N, gates, out)


def kernel_from_W(F, r, W, w):
    """explicit r-point DFT circuit from auxiliary points W (fs_search output)."""
    rows = [tuple(F.pw(w, k * n) for n in range(r)) for k in range(r)]
    have = {i: tuple(1 if c == i else 0 for c in range(r)) for i in range(r)}
    names, gates, nxt = {}, [], r
    todo = [("w", tuple(x)) for x in W] + [("X%d" % k, rows[k]) for k in range(r)]
    outs = [None] * r
    while todo:
        prog = False
        for item in list(todo):
            name, t = item
            for (iu, u), (iv, v) in itertools.combinations(have.items(), 2):
                sol = next(((a, b) for a in range(1, F.q) for b in range(1, F.q)
                            if F.lin(a, u, b, v) == t), None)
                if sol:
                    gates.append((sol[0], iu, sol[1], iv)); have[nxt] = t
                    if name.startswith("X"):
                        outs[int(name[1:])] = nxt
                    nxt += 1; todo.remove(item); prog = True
                    break
        if not prog:
            raise RuntimeError("cannot rebuild kernel")
    return (r, gates, outs)


def five_point(F, w):
    """Lev's 15-gate and Winograd's 17-gate 5-point DFT circuits over GF(q) (Lev 1980,
    Sec. 5.4.3; checked over Q(omega) by verify_witnesses.py).  Nodes 0..4 are inputs."""
    M, A, N_ = F.mul, F.add, F.neg
    half = F.inv[2]
    pw = [F.pw(w, k) for k in range(5)]
    g = M[half][A[pw[1]][pw[4]]]                      # cos u
    h = M[half][A[pw[2]][pw[3]]]                      # cos 2u
    f = M[half][A[pw[1]][N_[pw[4]]]]                  # i sin u
    d = M[half][A[pw[2]][N_[pw[3]]]]                  # i sin 2u
    hi = F.inv[h]
    one, m1 = 1, N_[1]

    def build(spec):
        ids, gates = {"y%d" % i: i for i in range(5)}, []
        for name, a, u, b, v in spec:
            gates.append((a, ids[u], b, ids[v])); ids[name] = 4 + len(gates)
        return (5, gates, [ids["X%d" % k] for k in range(5)])
    common = [("S1", one, "y1", one, "y4"), ("S2", one, "y1", m1, "y4"),
              ("S3", one, "y3", one, "y2"), ("S4", one, "y3", m1, "y2")]
    lev = build(common + [
        ("E1", one, "y0", h, "S1"), ("E2", one, "E1", h, "S3"),
        ("X0", hi, "E2", A[one][N_[hi]], "y0"),
        ("S10", one, "E2", A[g][N_[h]], "S1"), ("S11", one, "E1", g, "S3"),
        ("S12", f, "S2", N_[d], "S4"), ("S13", d, "S2", f, "S4"),
        ("X1", one, "S10", one, "S12"), ("X4", one, "S10", m1, "S12"),
        ("X2", one, "S11", one, "S13"), ("X3", one, "S11", m1, "S13")])
    a = A[M[half][A[g][h]]][m1]
    b = M[half][A[g][N_[h]]]
    c, e = A[f][d], A[f][N_[d]]
    win = build(common + [
        ("S5", one, "S1", one, "S3"), ("S6", one, "S1", m1, "S3"), ("S7", one, "S2", one, "S4"),
        ("X0", one, "S5", one, "y0"), ("S9", one, "X0", a, "S5"),
        ("S10", one, "S9", b, "S6"), ("S11", one, "S9", N_[b], "S6"),
        ("S12", c, "S2", N_[d], "S7"), ("S13", d, "S7", e, "S4"),
        ("X1", one, "S10", one, "S12"), ("X4", one, "S10", m1, "S12"),
        ("X2", one, "S11", one, "S13"), ("X3", one, "S11", m1, "S13")])
    return {"Lev 15": lev, "Winograd 17": win}


def values(F, C):
    n, gates, _ = C
    val = [tuple(1 if c == i else 0 for c in range(n)) for i in range(n)]
    for a, u, b, v in gates:
        val.append(F.lin(a, val[u], b, val[v]))
    return val


def check_dft(F, C, w):
    n = C[0]
    val = values(F, C)
    return all(val[C[2][k]] == tuple(F.pw(w, k * m) for m in range(n)) for k in range(n))


# ---------------------------------------------------------------- windows
def windows(C, smax):
    n, gates, outs = C
    G = len(gates)
    nodes = n + G
    preds = {n + g: (u, v) for g, (_, u, _, v) in enumerate(gates)}
    succ = {i: set() for i in range(nodes)}
    for x, (u, v) in preds.items():
        succ[u].add(x); succ[v].add(x)
    desc = {}
    for x in reversed(range(nodes)):
        d = 0
        for y in succ[x]:
            d |= (1 << y) | desc[y]
        desc[x] = d
    gate_ids = list(range(n, nodes))
    seen, frontier = set(), [1 << g for g in gate_ids]
    seen.update(frontier)
    while frontier:
        new = []
        for S in frontier:
            yield S
            if bin(S).count("1") >= smax:
                continue
            members = [x for x in gate_ids if S >> x & 1]
            nb = set()
            for x in members:
                nb.update(y for y in preds[x] if y >= n)
                nb.update(succ[x])
            for y in nb:
                T = S | (1 << y)
                if T not in seen:
                    seen.add(T); new.append(T)
        frontier = new


def is_convex(C, S, desc_cache={}):
    n, gates, _ = C
    key = id(C)
    if key not in desc_cache:
        nodes = n + len(gates)
        succ = {i: set() for i in range(nodes)}
        for g, (_, u, _, v) in enumerate(gates):
            succ[u].add(n + g); succ[v].add(n + g)
        desc = {}
        for x in reversed(range(nodes)):
            d = 0
            for y in succ[x]:
                d |= (1 << y) | desc[y]
            desc[x] = d
        desc_cache[key] = desc
    desc = desc_cache[key]
    reach_out = 0
    for x in range(n + len(gates)):
        if S >> x & 1:
            reach_out |= desc[x]
    outside = reach_out & ~S
    for x in range(n + len(gates)):
        if outside >> x & 1 and desc[x] & S:
            return False
    return True


def window_problem(F, C, val, S):
    n, gates, outs = C
    used_outside = set(outs)
    I, O = [], []
    for g, (_, u, _, v) in enumerate(gates):
        x = n + g
        if S >> x & 1:
            for y in (u, v):
                if not S >> y & 1 and y not in I:
                    I.append(y)
        else:
            for y in (u, v):
                if S >> y & 1:
                    used_outside.add(y)
    for x in range(n, n + len(gates)):
        if S >> x & 1 and x in used_outside:
            O.append(x)
    return I, O


def run_fs(d, q, targets, wmax, maxnodes):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write("%d\n" % len(targets))
        for t in targets:
            f.write(" ".join(map(str, t)) + "\n")
        path = f.name
    try:
        res = subprocess.run([os.path.join(HERE, "fs_search"), str(d), str(q), str(wmax), str(wmax),
                              "file:" + path, "maxnodes", str(maxnodes)],
                             capture_output=True, text=True)
        out = res.stdout
    finally:
        os.unlink(path)
    # strict: only a successful run with a recognized terminal line counts as a result
    if res.returncode != 0:
        raise RuntimeError(f"fs_search failed (exit {res.returncode}): {res.stderr.strip()[:500]}")
    if "FOUND" in out:
        return "FOUND", out
    if "ABORT" in out:
        return "ABORT", out
    if re.search(r"\|W\|<=\d+ \(L<=\d+\): none\s+nodes=\d+", out):
        return "none", out
    raise RuntimeError(f"fs_search gave no recognized result: {out.strip()[-500:]!r} {res.stderr.strip()[:300]!r}")


def main():
    q = int(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    DMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    SMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    MAXN = int(float(sys.argv[5])) if len(sys.argv) > 5 else 20_000_000
    WORKERS = int(sys.argv[6]) if len(sys.argv) > 6 else 12
    F = GF(q)
    w = F.root(N)
    if N == 5:
        circuits = five_point(F, w)
        print(f"GF({q}) N=5 omega={w}; starting circuits {list(circuits)}; DMAX={DMAX} SMAX={SMAX}", flush=True)
        return search(F, q, w, circuits, DMAX, SMAX, MAXN, WORKERS)
    w4 = F.pw(w, N // 4)
    K2 = kernel2(F)
    K4 = [compose(F, 2, K2, 2, K2, w4)]
    # extra 4-point kernels: all distinct 8-gate circuits from fs_search on this field
    out = subprocess.run([os.path.join(HERE, "fs_search"), "4", str(q), "4", "4", "dft", "nosym", "all", "12"],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.strip().startswith("W:"):
            Wpts = [tuple(map(int, t.strip("()").split(","))) for t in line.split(":")[1].split()]
            try:
                k = kernel_from_W(F, 4, Wpts, w4)
                if check_dft(F, k, w4):
                    K4.append(k)
            except RuntimeError:
                pass
    circuits = {}
    for i, k4 in enumerate(K4):
        if N == 4:
            circuits["K4#%d" % i] = k4
        else:
            circuits["DIT 4x2 (K4#%d)" % i] = compose(F, 4, k4, 2, K2, w)
            circuits["DIF 2x4 (K4#%d)" % i] = compose(F, 2, K2, 4, k4, w)
    print(f"GF({q}) N={N} omega={w}; {len(K4)} four-point kernels; DMAX={DMAX} SMAX={SMAX}", flush=True)
    search(F, q, w, circuits, DMAX, SMAX, MAXN, WORKERS)


def search(F, q, w, circuits, DMAX, SMAX, MAXN, WORKERS):
    cache, stats = {}, {"windows": 0, "convex": 0, "problems": 0, "aborted": 0, "found": 0}
    for name, C in circuits.items():
        assert check_dft(F, C, w), name
        val = values(F, C)
        nwin = 0
        for S in windows(C, SMAX):
            stats["windows"] += 1
            if not is_convex(C, S):
                continue
            stats["convex"] += 1
            I, O = window_problem(F, C, val, S)
            size = bin(S).count("1")
            if len(I) > DMAX or F.rank([val[i] for i in I]) < len(I):
                continue
            coords = [F.solve([val[i] for i in I], val[o]) for o in O]
            tg = sorted({F.normalize(tuple(c)) for c in coords})
            tg = [t for t in tg if sum(1 for x in t if x) > 1]
            wmax = size - 1 - len(tg)
            if wmax < 0:
                continue
            key = (len(I), tuple(tg), wmax)
            if key in cache:
                continue
            stats["problems"] += 1; nwin += 1
            cache[key] = (name, size)
        print(f"  {name}: {nwin} new window problems", flush=True)
    dump = os.environ.get("DUMP")
    if dump:  # write the window problems as fs_search input files instead of solving them
        os.makedirs(dump, exist_ok=True)
        with open(os.path.join(dump, "manifest.txt"), "w") as man:
            for i, ((d, tg, wmax), (name, size)) in enumerate(cache.items()):
                fn = "prob%04d.txt" % i
                with open(os.path.join(dump, fn), "w") as f:
                    f.write("%d\n" % len(tg))
                    for t in tg:
                        f.write(" ".join(map(str, t)) + "\n")
                man.write(f"{fn} q={q} d={d} wmax={wmax} start={name.replace(' ', '_')} size={size}\n")
        print(f"dumped {len(cache)} problems to {dump}", flush=True)
        return
    from concurrent.futures import ThreadPoolExecutor
    def job(item):
        (d, tg, wmax), (name, size) = item
        return item, run_fs(d, q, list(tg), wmax, MAXN)
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for ((d, tg, wmax), (name, size)), (res, txt) in ex.map(job, cache.items()):
            if res == "ABORT":
                stats["aborted"] += 1
                print(f"  abort: d={d} wmax={wmax} ({name}, {size} gates)", flush=True)
            if res == "none":
                stats["none"] = stats.get("none", 0) + 1
            if res == "FOUND":
                stats["found"] += 1
                print(f"*** IMPROVEMENT in {name}: window {size} gates, d={d}, targets={list(tg)}\n{txt}", flush=True)
    print("summary:", stats, flush=True)


if __name__ == "__main__":
    main()
