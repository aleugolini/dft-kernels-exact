#!/usr/bin/env python3
"""Independent cross-check of the smallest exhaustive searches (no code shared with
groebner/exact_small.py, no Singular).

A circuit with G gates is enumerated in the most literal way: gate g = a_g*u + b_g*v,
with u < v any two earlier nodes and BOTH coefficients unknown (no normalization), and
each gate is either one of the targets (exactly equal to that row) or an auxiliary
value.  No ordering rules, no support pruning, no symmetry: every sequence is visited.
After each target, sympy computes a Groebner basis over Q of the accumulated equations
together with Phi_N(w) (w is an unknown; since Phi_N is irreducible, the system is
inconsistent over Q(w)-bar iff the basis is [1]).  A leaf has all targets placed and
every auxiliary gate used.  No leaf proves that no circuit with <= G gates exists over
any field of characteristic 0 containing a primitive N-th root of unity.

usage: python3 independent_check.py KERNEL G [--split k P]
       KERNEL G in: dft3 4 | dft3 5 | wht4 6 | dft4 7
       --split k P: keep only the first gates whose parent pair has index = k mod P
"""
import sys
import time

import sympy as sp


def kernel(name):
    w = sp.Symbol("w")
    if name.startswith("dft"):
        N = int(name[3:])
        rows = [[w ** ((k * n) % N) for n in range(N)] for k in range(N)]
        return N, rows, [sp.cyclotomic_poly(N, w)], w
    if name == "wht4":
        rows = [[(-1) ** bin(k & n).count("1") for n in range(4)] for k in range(4)]
        return 4, rows, [], None


def main():
    name, G = sys.argv[1], int(sys.argv[2])
    split = (0, 1)
    if "--split" in sys.argv:
        i = sys.argv.index("--split")
        split = (int(sys.argv[i + 1]), int(sys.argv[i + 2]))
    N, rows, base, w = kernel(name)
    K = len(rows)
    a = sp.symbols(f"a0:{G}")
    b = sp.symbols(f"b0:{G}")
    gens = list(a) + list(b) + ([w] if w is not None else [])
    stats = dict(nodes=0, gb=0, leaves=0)
    forms = [[sp.Integer(1 if c == j else 0) for c in range(N)] for j in range(N)]
    parents = []      # parent pairs of the gates
    role = []         # target index or None for each gate

    def consistent(eqs):
        stats["gb"] += 1
        gb = sp.groebner(eqs, *gens, order="grevlex", domain=sp.QQ)
        return list(gb.exprs) != [1]

    def dfs(eqs, done):
        stats["nodes"] += 1
        g = len(parents)
        if len(done) == K:
            used = {p for pr in parents for p in pr}
            if all(r is not None or N + i in used for i, r in enumerate(role)):
                stats["leaves"] += 1
            return
        if g + (K - len(done)) > G:
            return
        pos = N + g
        pairs = [(u, v) for u in range(pos) for v in range(u + 1, pos)]
        for idx, (u, v) in enumerate(pairs):
            if g == 0 and idx % split[1] != split[0]:
                continue
            form = [sp.expand(a[g] * forms[u][c] + b[g] * forms[v][c]) for c in range(N)]
            parents.append((u, v))
            forms.append(form)
            for k in range(K):  # gate g is target k
                if k in done:
                    continue
                new = eqs + [form[c] - rows[k][c] for c in range(N)]
                if consistent(new):
                    role.append(k)
                    dfs(new, done | {k})
                    role.pop()
            if g + 1 + (K - len(done)) <= G:  # gate g is auxiliary
                role.append(None)
                dfs(eqs, done)
                role.pop()
            forms.pop()
            parents.pop()

    t0 = time.time()
    dfs(list(base), frozenset())
    stats["sec"] = round(time.time() - t0, 1)
    print(f"{name}, at most {G} gates, unnormalized, no ordering rules, sympy, split={split}: {stats}")
    print("RESULT:", "circuits exist" if stats["leaves"] else f"no circuit -> L > {G}")


if __name__ == "__main__":
    main()
