#!/usr/bin/env python3
"""Step 2: exact free-scalar gate count of small kernels over C (or any algebraically
closed field), by exhaustive enumeration of circuit wirings and Groebner bases.

Model (D): every gate is z = a*u + b*v with free scalars.  A *minimal* circuit has no
gate that is a multiple of an earlier node and no unused non-output gate, and computes
each target (a row of the matrix, up to scale) exactly once.  So it is a sequence of
nodes after the inputs e_0..e_{N-1}:
  W node (auxiliary):  w = p + t*q          (p, q earlier nodes, t a new unknown)
  T node (target k):   T_k = a_k*p + b_k*q   (a_k, b_k unknowns, T_k the exact row)
Give each node the key ('T', k) or ('W', parents); every minimal circuit has a
linear extension that is lexicographically smallest in these keys, and it satisfies:
(i) a target that follows a W node uses it, (ii) adjacent independent targets appear
in increasing k, (iii) adjacent independent W nodes have non-decreasing parent pairs.
With --sym2 (N prime) the first gate is fixed to e_0 + t e_1: the affine maps
n -> an+b permute the DFT rows up to scale and act 2-transitively on the inputs, and
the smallest-key first node is an input-input W with the smallest pair.  The search runs
over all such sequences with at most G gates, with two exact prunings:
  * support: supp(T_k) must lie in supp(p) U supp(q) (generic supports are upper bounds);
  * algebra: after every target, std() of the accumulated equations over Q(omega)
    (Singular); if it is (1) no completion has a solution over the algebraic closure.
Singular is taken from $SINGULAR if set (e.g. a micromamba env on the servers).
A surviving leaf is a wiring whose system is consistent: a circuit with <= G gates.
No surviving leaf proves L(kernel) > G over every field of characteristic 0 that
contains omega (in particular over C).

usage: python3 exact_small.py KERNEL G [--char p --omega v] [--show K] [--sym2]
                               [--split k K] [--split-depth d]   (branches at gate d=3, k mod K)
  KERNEL: dftN (N = 2..8) or wht4
  --char p --omega v : run over F_p-bar with omega = v in F_p (controls)
  --no-rules / --no-support : disable R1-R3 / the support pruning (cross-checks)
  --over-z : Groebner bases over Z with omega an unknown and Phi_N(omega) in the ideal; every
             pruned prefix records its integer obstruction d (I cap Z = (d)); a circuit with
             <= G gates over a field of characteristic p exists only if p divides some d
             (or a leaf exists)
  --galois : (dftN, N prime, char 0) the automorphisms omega -> omega^j of Q(omega) map a
             circuit to one with row labels k -> jk (mod N), so the first target with a
             nonzero row index can be required to be row 1 (orbit of size N-1)
  --cert FILE : (characteristic 0 over Q(omega), or --char p) for every branch discarded by the
             Groebner test, write a certificate: the polynomials h_i with sum_i h_i f_i = 1,
             where f_i are the equations of the prefix in order (Singular's lift), one JSON line
             per branch (gzip if FILE ends in .gz); checked by src/check_certificates.py
Every PROGRESS_EVERY seconds (default 600, 0 = off) a line "progress: ..." goes to stderr,
with the estimated fraction of this part done (from the choices made at the top 4 levels).
"""
import gzip
import hashlib
import json
import os
import resource
import subprocess
import sys
import time

CYCLO = {2: "w+1", 3: "w2+w+1", 4: "w2+1", 5: "w4+w3+w2+w+1", 6: "w2-w+1",
         7: "w6+w5+w4+w3+w2+w+1", 8: "w4+1"}


class Singular:
    """Persistent Singular process.  Singular's memory grows by about 2 KB per command in
    this workload, so the process is restarted every `restart_every` commands and the
    current search state (self.state: the commands that define it) is replayed."""

    def __init__(self, restart_every=50_000):
        self.restart_every = restart_every
        self.state = []
        self.calls = 0
        self.restarts = 0
        self.p = None
        self._start()

    def _start(self):
        if self.p is not None:
            self.p.kill()
            self.p.wait()
            self.restarts += 1
        self.p = subprocess.Popen([os.environ.get("SINGULAR", "Singular"), "-q", "--no-warn"],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, text=True, bufsize=1)
        self.n = 0
        for cmd in self.state:  # replay: same check as push(), state commands print nothing
            out = self._exec(cmd)
            if out:
                raise RuntimeError("unexpected Singular output on replay:\n" + "\n".join(out) + "\ncommand: " + cmd[:500])

    def _exec(self, cmd):
        self.n += 1
        self.p.stdin.write(cmd + '\nprint("@@END");\n')
        self.p.stdin.flush()
        out = []
        while True:
            line = self.p.stdout.readline()
            if not line:
                raise RuntimeError("Singular died:\n" + "\n".join(out))
            if line.strip() == "@@END":
                bad = [l for l in out if l.lstrip().startswith("?") or "error" in l.lower()]
                if bad:  # an error of Singular must never be read as a result (strict protocol)
                    raise RuntimeError("Singular error:\n" + "\n".join(out) + "\ncommand: " + cmd[:500])
                return out
            out.append(line.rstrip())

    def run(self, cmd):
        """Run a command that does not change the search state."""
        if self.n >= self.restart_every:
            self._start()
        self.calls += 1
        return self._exec(cmd)

    def push(self, cmd):
        """Run a command that becomes part of the search state (it must print nothing)."""
        out = self.run(cmd)
        if out:
            raise RuntimeError("unexpected Singular output:\n" + "\n".join(out) + "\ncommand: " + cmd[:500])
        self.state.append(cmd)
        return out

    def close(self):
        """Stop the process and wait for it (so that its CPU time is accounted)."""
        if self.p is not None:
            self.p.stdin.close()
            self.p.wait()
            self.p = None

    def pop(self, k=1):
        del self.state[len(self.state) - k:]


CYCLO_Z = {3: "w^2+w+1", 4: "w^2+1", 5: "w^4+w^3+w^2+w+1", 6: "w^2-w+1",
           7: "w^6+w^5+w^4+w^3+w^2+w+1", 8: "w^4+1"}


def setup(kernel, G, char, omega, overz=False):
    if kernel.startswith("dft"):
        N = int(kernel[3:])
        w = "w" if char == 0 else str(omega)
        targets = [[f"({w})^{(k * n) % N}" for n in range(N)] for k in range(N)]
    elif kernel == "wht4":
        N = 4
        targets = [[str((-1) ** bin(k & n).count("1")) for n in range(N)] for k in range(N)]
    else:
        raise SystemExit("unknown kernel")
    K = len(targets)
    nW = G - K
    vars_ = [f"t{j}" for j in range(1, nW + 1)] + [f"a{k}" for k in range(K)] + \
        [f"b{k}" for k in range(K)]
    if overz:  # integers, omega an unknown, Phi_N(omega) added to every ideal
        ring = f"ring R = (integer), (w,{','.join(vars_)}), dp;"
    elif char == 0 and kernel.startswith("dft") and N > 2:
        ring = f"ring R = (0,w), ({','.join(vars_)}), dp; minpoly = {CYCLO[N]};"
    else:
        ring = f"ring R = ({char}), ({','.join(vars_)}), dp;"
        targets = [[t.replace("(w)", "(-1)") for t in row] for row in targets]
    return N, targets, nW, ring


def search(kernel, G, char=0, omega=None, show=3, sym2=False, split=(0, 1), split_depth=3,
           rules=True, support=True, galois=False, overz=False, cert=None, manifest=None, only_path=None):
    """only_path: optional list of moves (kind, k, u, v); the search then follows only this sequence,
    with all the rules and tests active (positive control: a known circuit must give a leaf)."""
    manifest = {} if manifest is None else manifest
    N, targets, nWmax, ring = setup(kernel, G, char, omega, overz)
    K = len(targets)
    S = Singular(int(os.environ.get("RESTART_EVERY", 50_000)))
    ver = S.run('print(string(system("version")));')
    manifest.update(singular_version=ver[0].strip() if ver else "?", restart_every=S.restart_every,
                    singular=os.environ.get("SINGULAR", "Singular"))
    if overz:
        S.push(ring + f" option(redSB); ideal G0 = std(ideal({CYCLO_Z[N]})); int ii;")
    else:
        S.push(ring + " option(redSB); ideal G0 = 0;")
    for j in range(N):
        S.push("; ".join(f"poly n{j}_{c} = {1 if c == j else 0}" for c in range(N)) + ";")
    for k in range(K):
        S.push("; ".join(f"poly T{k}_{c} = {targets[k][c]}" for c in range(N)) + ";")
    full = [frozenset(c for c in range(N) if targets[k][c] not in ("0",)) for k in range(K)]

    nodes = [dict(kind="E", k=j, par=(), supp=frozenset([j])) for j in range(N)]
    stats = dict(dfs=0, gb=0, leaves=0)
    # per gate depth d (1 = first gate): moves generated, Groebner tests, tests giving (1)
    depth = {}
    branch = [-1]
    sols = []

    obstr = {}  # over Z: integer d of every pruned prefix -> count
    eqstack = []  # equations of the targets on the current path, in order (for --cert)
    certf = (gzip.open(cert, "wt") if cert.endswith(".gz") else open(cert, "w")) if cert else None

    def write_cert(gens, move):
        out = S.run(f"ideal I_c = {','.join(gens)}; matrix T_c = lift(I_c, ideal(1));"
                    f" for (int i_c = 1; i_c <= nrows(T_c); i_c++) {{ print(string(T_c[i_c, 1])); }}"
                    f" kill I_c; kill T_c; kill i_c;")
        if len(out) != len(gens):
            raise RuntimeError(f"lift returned {len(out)} cofactors for {len(gens)} equations")
        prefix = [[nd["kind"], nd["k"], *nd["par"]] for nd in nodes[N:]] + [list(move)]
        certf.write(json.dumps({"p": prefix, "h": out}, separators=(",", ":")) + "\n")

    def add_and_check(stdcmd, dst):
        stats["gb"] += 1
        if overz:  # the constant of the strong Groebner basis over Z generates I cap Z
            out = S.run(stdcmd + f" poly c = 0; for (ii = 1; ii <= size(G{dst}); ii++)"
                        f" {{ if (deg(G{dst}[ii]) == 0 && G{dst}[ii] != 0) {{ c = G{dst}[ii]; }} }} print(c);")
            if len(out) != 1 or not out[0].strip().lstrip("-").isdigit():
                raise RuntimeError(f"unexpected Singular output for {stdcmd[:200]}: {out}")
            d = abs(int(out[0].strip()))
            if d:
                obstr[d] = obstr.get(d, 0) + 1
            return d == 0
        out = S.run(stdcmd + f" int f = (deg(G{dst}[1]) == 0 && G{dst}[1] != 0); print(f); kill f;")
        if out not in (["0"], ["1"]):  # only an explicit 0/1 answer is accepted
            raise RuntimeError(f"unexpected Singular output for {stdcmd[:200]}: {out}")
        return out == ["0"]

    path = []  # (index, count) of the choices made at the top levels, for progress
    every = float(os.environ.get("PROGRESS_EVERY", 600))
    clock = dict(t0=time.time(), last=time.time())

    def progress():
        frac, w = 0.0, 1.0
        for i, n in path[:4]:
            frac += w * i / n
            w /= n
        now = time.time()
        print(f"progress: {now - clock['t0']:.0f} s, ~{100 * frac:.1f}% of this part, "
              f"gb={stats['gb']}, path={' '.join(f'{i + 1}/{n}' for i, n in path[:4])}",
              file=sys.stderr, flush=True)
        clock["last"] = now

    def dfs(gidx, nW, done):
        stats["dfs"] += 1
        if every > 0 and time.time() - clock["last"] >= every:
            progress()
        pos = len(nodes)
        nG = pos - N
        if len(done) == K:
            used = {p for nd in nodes for p in nd["par"]}
            if all(i in used for i, nd in enumerate(nodes) if nd["kind"] == "W"):
                stats["leaves"] += 1
                if len(sols) < show:
                    gb = S.run(f"print(G{gidx});")
                    sols.append(([(nd["kind"], nd["k"], nd["par"]) for nd in nodes[N:]], gb))
            return
        if nG + (K - len(done)) > G:
            return
        last = nodes[-1]

        def keep():  # work split over processes, counted on both kinds of branches
            if nG + 1 == split_depth:
                branch[0] += 1
                return branch[0] % split[1] == split[0]
            return True

        moves = []
        # targets
        for k in range(K):
            if k in done:
                continue
            if galois and k not in (0, 1) and not (done - {0}):
                continue  # Galois: the first target with a nonzero row index is row 1
            for u in range(pos):
                for v in range(u + 1, pos):
                    if rules and last["kind"] == "W" and pos - 1 not in (u, v):
                        continue  # R1
                    if rules and last["kind"] == "T" and pos - 1 not in (u, v) and k < last["k"]:
                        continue  # R2
                    if support and not full[k] <= nodes[u]["supp"] | nodes[v]["supp"]:
                        continue
                    if keep():
                        moves.append(("T", k, u, v))
        # auxiliary nodes
        if nW < nWmax and nG + 1 + (K - len(done)) <= G:
            for u in range(pos):
                for v in range(u + 1, pos):  # p + t q, t free: the line minus q (a copy)
                    if nG == 0 and sym2 and (u, v) != (0, 1):
                        continue
                    if rules and last["kind"] == "W" and pos - 1 not in (u, v) and (u, v) < last["par"]:
                        continue  # R3
                    if keep():
                        moves.append(("W", None, u, v))

        if only_path is not None:  # positive control: keep only the next move of the given sequence
            nxt = tuple(only_path[nG]) if nG < len(only_path) else None
            moves = [m for m in moves if m == nxt]
        dd = depth.setdefault(nG + 1, [0, 0, 0])
        dd[0] += len(moves)
        path.append(None)
        for i, (kind, k, u, v) in enumerate(moves):
            path[-1] = (i, len(moves))
            if kind == "T":
                eqs = ",".join(f"a{k}*n{u}_{c}+b{k}*n{v}_{c}-T{k}_{c}" for c in range(N))
                stdcmd = f"ideal G{gidx + 1} = std(G{gidx}, ideal({eqs}));"
                dd[1] += 1
                if not add_and_check(stdcmd, gidx + 1):
                    dd[2] += 1
                    if certf:
                        write_cert(eqstack + [f"a{k}*n{u}_{c}+b{k}*n{v}_{c}-T{k}_{c}" for c in range(N)],
                                   ("T", k, u, v))
                    continue
                S.state.append(stdcmd)
                eqstack.extend(f"a{k}*n{u}_{c}+b{k}*n{v}_{c}-T{k}_{c}" for c in range(N))
                S.push("; ".join(f"poly n{pos}_{c} = T{k}_{c}" for c in range(N)) + ";")
                nodes.append(dict(kind="T", k=k, par=(u, v), supp=full[k]))
                dfs(gidx + 1, nW, done | {k})
                nodes.pop()
                S.pop(2)
                del eqstack[len(eqstack) - N:]
            else:
                t = f"t{nW + 1}"
                S.push("; ".join(f"poly n{pos}_{c} = n{u}_{c}+{t}*n{v}_{c}" for c in range(N)) + ";")
                nodes.append(dict(kind="W", k=None, par=(u, v),
                                  supp=nodes[u]["supp"] | nodes[v]["supp"]))
                dfs(gidx, nW + 1, done)
                nodes.pop()
                S.pop()
        path.pop()

    t0 = time.time()
    dfs(0, 0, frozenset())
    stats["sec"] = round(time.time() - t0, 1)
    S.close()  # CPU time of this process and of all the Singular processes (user + system)
    ru_s, ru_c = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    stats["cpu"] = round(ru_s.ru_utime + ru_s.ru_stime + ru_c.ru_utime + ru_c.ru_stime, 1)
    stats["restarts"] = S.restarts
    stats["depth"] = {d: tuple(v) for d, v in sorted(depth.items())}
    stats["obstr"] = dict(sorted(obstr.items()))
    if certf:
        certf.close()
    return stats, sols


def main():
    kernel, G = sys.argv[1], int(sys.argv[2])
    char = int(sys.argv[sys.argv.index("--char") + 1]) if "--char" in sys.argv else 0
    omega = int(sys.argv[sys.argv.index("--omega") + 1]) if "--omega" in sys.argv else None
    show = int(sys.argv[sys.argv.index("--show") + 1]) if "--show" in sys.argv else 3
    sym2 = "--sym2" in sys.argv
    split = (0, 1)
    if "--split" in sys.argv:
        i = sys.argv.index("--split")
        split = (int(sys.argv[i + 1]), int(sys.argv[i + 2]))
    sd = int(sys.argv[sys.argv.index("--split-depth") + 1]) if "--split-depth" in sys.argv else 3
    rules = "--no-rules" not in sys.argv
    support = "--no-support" not in sys.argv
    galois = "--galois" in sys.argv
    if galois:
        n = int(kernel[3:]) if kernel.startswith("dft") else 0
        if char != 0 or n < 3 or any(n % d == 0 for d in range(2, n)):
            raise SystemExit("--galois needs dftN with N prime >= 3 over characteristic 0")
    overz = "--over-z" in sys.argv
    cert = sys.argv[sys.argv.index("--cert") + 1] if "--cert" in sys.argv else None
    if cert and (overz or galois or split != (0, 1)):
        raise SystemExit("--cert is for unsplit runs over Q(omega) or F_p (no --over-z, --galois)")
    if overz and (char != 0 or not kernel.startswith("dft") or galois):
        raise SystemExit("--over-z needs dftN in characteristic 0, without --galois")
    src = os.path.abspath(__file__)
    manifest = dict(source=os.path.basename(src), source_sha256=hashlib.sha256(open(src, "rb").read()).hexdigest(),
                    argv=sys.argv[1:], split_depth=sd, python=sys.version.split()[0])
    stats, sols = search(kernel, G, char, omega, show, sym2, split, sd, rules, support, galois, overz, cert,
                         manifest)
    print("manifest:", json.dumps(manifest, sort_keys=True))
    per_depth = stats.pop("depth")
    obstr = stats.pop("obstr")
    field = ("Z (Phi_N in the ideal)" if overz else "Q(omega)-bar") if char == 0 else f"F_{char}-bar (omega={omega})"
    extra = "" if rules and support else f", rules={rules}, support={support}"
    extra += ", galois=True" if galois else ""
    print(f"{kernel}, at most {G} gates, over {field}, sym2={sym2}, split={split}{extra}: {stats}")
    print("per depth (gate: moves, Groebner tests, rejected):",
          "; ".join(f"{d}: {m}, {g}, {r}" for d, (m, g, r) in per_depth.items()))
    print("RESULT:", "circuits exist" if stats["leaves"] else f"no circuit -> L > {G}")
    if overz:
        primes = sorted({p for d in obstr for p in range(2, d + 1)
                         if d % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))})
        print("obstructions d (value: pruned prefixes):", obstr)
        print("candidate characteristics (prime divisors of the d's):", primes)
    for wiring, gb in sols:
        print("wiring:", wiring)
        print("  GB:", " ".join(gb))


if __name__ == "__main__":
    main()
