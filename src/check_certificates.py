#!/usr/bin/env python3
"""Independent check of the negative results of groebner/exact_small.py with certificates.

The program shares no code with the search and does not use Singular or any computer
algebra system. It enumerates, with its own code, the sequences of auxiliary and target
nodes of Section III of the letter (ordering rules R1-R3, support test, budget of G gates,
and, with --sym2, first node x0 + t x1), and builds the polynomial equations of every
target move itself. A target move is discarded only if the certificate file contains, for
exactly that prefix, polynomials h_i such that

    sum_i h_i f_i = 1   in Q(omega)[t, a, b]  (or F_p[t, a, b] with --char p),

where f_i are the equations of the prefix (the N coordinates of a_k u + b_k v - f_k, for
the targets in order); the identity is verified exactly, with rational arithmetic and
reduction modulo the cyclotomic polynomial. Any other target move is followed. The bound
L > G is confirmed only if every sequence is discarded, i.e., no leaf is reached; a missing
certificate can only make the check fail, never pass. The certificates are written by
`exact_small.py ... --cert FILE`, one JSON line per discarded branch.

usage: python3 check_certificates.py KERNEL G CERTFILE [--sym2] [--char p --omega v] [--stream]
       KERNEL: dftN (N = 3..8) or wht4
       CERTFILE: a file, or a pattern such as 'dir/cert.part*.jsonl.gz' (files read in sorted order)
       --stream: read the certificates one at a time, in the order of the enumeration, instead of
                 loading them all (constant memory, for large files); a certificate out of order
                 is not found, which can only make the check fail
"""
import glob
import gzip
import json
import re
import sys
from fractions import Fraction

CYCLO = {3: [1, 1, 1], 4: [1, 0, 1], 5: [1, 1, 1, 1, 1], 6: [1, -1, 1],
         7: [1] * 7, 8: [1, 0, 0, 0, 1]}  # coefficients of Phi_N, constant term first


class Ring:
    """Polynomials as dicts {exponent tuple: coefficient}; variable 0 is omega."""

    def __init__(self, names, char):
        self.names, self.char = names, char
        self.index = {n: i for i, n in enumerate(names)}
        self.n = len(names)

    def norm(self, c):
        return c % self.char if self.char else c

    def const(self, c):
        c = self.norm(c)
        return {(0,) * self.n: c} if c else {}

    def var(self, name, e=1):
        m = [0] * self.n
        m[self.index[name]] = e
        return {tuple(m): 1}

    def add(self, p, q, s=1):
        r = dict(p)
        for m, c in q.items():
            v = self.norm(r.get(m, 0) + s * c)
            if v:
                r[m] = v
            else:
                r.pop(m, None)
        return r

    def mul(self, p, q):
        r = {}
        for m1, c1 in p.items():
            for m2, c2 in q.items():
                m = tuple(a + b for a, b in zip(m1, m2))
                v = self.norm(r.get(m, 0) + c1 * c2)
                if v:
                    r[m] = v
                else:
                    r.pop(m, None)
        return r

    def pow(self, p, e):
        r = self.const(1)
        for _ in range(e):
            r = self.mul(r, p)
        return r

    def reduce_cyclo(self, p, phi):
        """Remainder modulo Phi_N(omega) (monic), as a polynomial in omega."""
        d = len(phi) - 1
        p = dict(p)
        while True:
            high = [m for m in p if m[0] >= d]
            if not high:
                return p
            m = max(high, key=lambda x: x[0])
            c = p.pop(m)
            for j in range(d):  # omega^d = -(phi_0 + ... + phi_{d-1} omega^{d-1})
                if phi[j]:
                    mm = (m[0] - d + j,) + m[1:]
                    v = self.norm(p.get(mm, 0) - c * phi[j])
                    if v:
                        p[mm] = v
                    else:
                        p.pop(mm, None)


TOKEN = re.compile(r"\s*(?:(\d+)|([A-Za-z]\w*)|(.))")


def parse(ring, s):
    """Parse a polynomial printed by Singular: integers, names, + - * / ^ and parentheses."""
    toks = [(a, b, c) for a, b, c in TOKEN.findall(s) if a or b or c.strip()]
    pos = [0]

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else (None, None, None)

    def take():
        pos[0] += 1
        return toks[pos[0] - 1]

    def atom():
        num, name, op = take()
        if num:
            return ring.const(int(num))
        if name:
            if name not in ring.index:
                raise ValueError(f"unknown variable {name!r}")
            return ring.var(name)
        if op == "(":
            e = expr()
            if take()[2] != ")":
                raise ValueError("missing )")
            return e
        if op == "-":
            return ring.mul(ring.const(-1), factor())
        raise ValueError(f"unexpected {op!r} in {s!r}")

    def factor():
        a = atom()
        if peek()[2] == "^":
            take()
            a = ring.pow(a, int(take()[0]))
        return a

    def term():
        a = factor()
        while peek()[2] in ("*", "/"):
            op = take()[2]
            b = factor()
            if op == "*":
                a = ring.mul(a, b)
            else:  # division by an integer constant only
                if len(b) != 1 or any(next(iter(b))):
                    raise ValueError("division by a non-constant")
                c = next(iter(b.values()))
                inv = pow(c, ring.char - 2, ring.char) if ring.char else Fraction(1, c)
                a = ring.mul(a, ring.const(inv))
        return a

    def expr():
        a = term() if peek()[2] != "-" else ring.mul(ring.const(-1), (take(), term())[1])
        while peek()[2] in ("+", "-"):
            op = take()[2]
            a = ring.add(a, term(), 1 if op == "+" else -1)
        return a

    e = expr()
    if pos[0] != len(toks):
        raise ValueError(f"trailing input in {s!r}")
    return e


def setup(kernel, G, char, omega):
    if kernel.startswith("dft"):
        N = int(kernel[3:])
        rows = [[(k * n) % N for n in range(N)] for k in range(N)]  # exponents of omega
    elif kernel == "wht4":
        N, rows = 4, None
    else:
        raise SystemExit("unknown kernel")
    K = N
    nW = G - K
    names = ["w"] + [f"t{j}" for j in range(1, nW + 1)] + [f"a{k}" for k in range(K)] + \
        [f"b{k}" for k in range(K)]
    ring = Ring(names, char)
    if rows is None:
        targets = [[ring.const((-1) ** bin(k & n).count("1")) for n in range(N)] for k in range(K)]
    elif char:
        targets = [[ring.const(pow(omega, e, char)) for e in row] for row in rows]
    else:
        targets = [[ring.var("w", e) if e else ring.const(1) for e in row] for row in rows]
    phi = CYCLO[N] if (kernel.startswith("dft") and not char) else None
    return N, K, nW, ring, targets, phi


def main():
    kernel, G, path = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    sym2 = "--sym2" in sys.argv
    char = int(sys.argv[sys.argv.index("--char") + 1]) if "--char" in sys.argv else 0
    omega = int(sys.argv[sys.argv.index("--omega") + 1]) if "--omega" in sys.argv else None
    N, K, nWmax, ring, targets, phi = setup(kernel, G, char, omega)
    supp = [frozenset(c for c in range(N) if targets[k][c]) for k in range(K)]

    stream = "--stream" in sys.argv
    files = sorted(glob.glob(path)) if any(ch in path for ch in "*?[") else [path]
    if not files:
        raise SystemExit(f"no certificate file matches {path}")

    def lines():
        for fn in files:
            with (gzip.open if fn.endswith(".gz") else open)(fn, "rt") as f:
                yield from f

    def decode(line):
        c = json.loads(line)
        return tuple(tuple(m) for m in c["p"]), c["h"]

    certs, nxt, nread = {}, [None], [0]
    src = lines()
    if not stream:
        for line in src:
            k, h = decode(line)
            certs[k] = h
        nread[0] = len(certs)

    def lookup(prefix):
        if not stream:
            return certs.get(prefix)
        if nxt[0] is None:
            line = next(src, None)
            if line is None:
                return None
            nxt[0] = decode(line)
            nread[0] += 1
        if nxt[0][0] == prefix:
            h, nxt[0] = nxt[0][1], None
            return h
        return None

    # nodes: (kind, k, (u, v), linear form = list of N polynomials, support)
    nodes = [("E", j, (), [ring.const(1 if c == j else 0) for c in range(N)], frozenset([j]))
             for j in range(N)]
    eqs = []  # equations of the targets on the current path
    stats = dict(nodes=0, target_moves=0, certified=0, leaves=0)
    used_certs = [0]

    def verify(prefix, gens):
        h = lookup(prefix)
        if h is None:
            return False
        if len(h) != len(gens):
            raise SystemExit(f"certificate for {prefix} has {len(h)} cofactors for {len(gens)} equations")
        total = ring.const(-1)
        for hs, g in zip(h, gens):
            if hs != "0":
                total = ring.add(total, ring.mul(parse(ring, hs), g))
        if phi:
            total = ring.reduce_cyclo(total, phi)
        if total:
            raise SystemExit(f"INVALID certificate for prefix {prefix}: sum h_i f_i - 1 != 0")
        used_certs[0] += 1
        return True

    def dfs(nW, done):
        stats["nodes"] += 1
        pos = len(nodes)
        nG = pos - N
        if len(done) == K:
            used = {p for nd in nodes for p in nd[2]}
            if all(i in used for i, nd in enumerate(nodes) if nd[0] == "W"):
                stats["leaves"] += 1
                print("LEAF (no certificate on this path):",
                      [(nd[0], nd[1], *nd[2]) for nd in nodes[N:]])
            return
        if nG + (K - len(done)) > G:
            return
        last = nodes[-1]
        prefix = tuple((nd[0], nd[1] if nd[0] == "T" else None, *nd[2]) for nd in nodes[N:])
        for k in range(K):  # target moves
            if k in done:
                continue
            for u in range(pos):
                for v in range(u + 1, pos):
                    if last[0] == "W" and pos - 1 not in (u, v):
                        continue  # R1
                    if last[0] == "T" and pos - 1 not in (u, v) and k < last[1]:
                        continue  # R2
                    if not supp[k] <= nodes[u][4] | nodes[v][4]:
                        continue  # support
                    stats["target_moves"] += 1
                    a, b = ring.var(f"a{k}"), ring.var(f"b{k}")
                    new = [ring.add(ring.add(ring.mul(a, nodes[u][3][c]), ring.mul(b, nodes[v][3][c])),
                                    targets[k][c], -1) for c in range(N)]
                    if verify(prefix + (("T", k, u, v),), eqs + new):
                        stats["certified"] += 1
                        continue
                    eqs.extend(new)
                    nodes.append(("T", k, (u, v), targets[k], supp[k]))
                    dfs(nW, done | {k})
                    nodes.pop()
                    del eqs[len(eqs) - N:]
        if nW < nWmax and nG + 1 + (K - len(done)) <= G:  # auxiliary moves
            t = ring.var(f"t{nW + 1}")
            for u in range(pos):
                for v in range(u + 1, pos):
                    if nG == 0 and sym2 and (u, v) != (0, 1):
                        continue
                    if last[0] == "W" and pos - 1 not in (u, v) and (u, v) < last[2]:
                        continue  # R3
                    form = [ring.add(nodes[u][3][c], ring.mul(t, nodes[v][3][c])) for c in range(N)]
                    nodes.append(("W", None, (u, v), form, nodes[u][4] | nodes[v][4]))
                    dfs(nW + 1, done)
                    nodes.pop()

    dfs(0, frozenset())
    if stream:  # count the certificates left in the files
        nread[0] += sum(1 for _ in src) + (nxt[0] is not None)
    unused = nread[0] - used_certs[0]
    field = f"F_{char} (omega={omega})" if char else ("Q(omega)" if phi else "Q")
    print(f"{kernel}, at most {G} gates, over {field}, sym2={sym2}: {stats}; certificates in file: "
          f"{nread[0]}, unused: {unused}")
    ok = not stats["leaves"]
    print("RESULT:", f"all branches certified -> L > {G}" if ok else "NOT CONFIRMED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
