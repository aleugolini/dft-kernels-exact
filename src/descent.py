#!/usr/bin/env python3
"""Iterated exact local improvement: start from a 24-gate F_N circuit, replace an
improvable window, rebuild and re-verify the whole circuit, repeat until no window
with <= DMAX inputs and <= SMAX gates can be shortened.

usage: python3 descent.py q start_circuit_name [N=8] [DMAX=6] [SMAX=10] [WORKERS=12]
       N=5: start_circuit_name is "Lev 15" or "Winograd 17"
"""
import itertools
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import window_search as ws


def base_circuits(F, N, w):
    if N == 5:
        return ws.five_point(F, w)
    w4 = F.pw(w, N // 4)
    K2 = ws.kernel2(F)
    K4 = [ws.compose(F, 2, K2, 2, K2, w4)]
    out = subprocess.run([os.path.join(ws.HERE, "fs_search"), "4", str(F.q), "4", "4", "dft", "nosym", "all", "12"],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.strip().startswith("W:"):
            Wp = [tuple(map(int, t.strip("()").split(","))) for t in line.split(":")[1].split()]
            k = ws.kernel_from_W(F, 4, Wp, w4)
            if ws.check_dft(F, k, w4):
                K4.append(k)
    C = {}
    for i, k4 in enumerate(K4):
        C["DIT 4x2 (K4#%d)" % i] = ws.compose(F, 4, k4, 2, K2, w)
        C["DIF 2x4 (K4#%d)" % i] = ws.compose(F, 2, K2, 4, k4, w)
    return C


def rebuild(F, N, rows, values_list):
    """explicit circuit (n, gates, outs) realising the given node values (any order)."""
    pending = list(dict.fromkeys(values_list))
    have = [tuple(1 if c == i else 0 for c in range(N)) for i in range(N)]
    gates = []
    while pending:
        prog = False
        for t in list(pending):
            for (iu, u), (iv, v) in itertools.combinations(list(enumerate(have)), 2):
                sol = next(((a, b) for a in range(1, F.q) for b in range(1, F.q) if F.lin(a, u, b, v) == t), None)
                if sol:
                    gates.append((sol[0], iu, sol[1], iv)); have.append(t); pending.remove(t); prog = True
                    break
        if not prog:
            return None
    return (N, gates, [have.index(r) for r in rows])


def check_random(F, C, rows):
    n, gates, outs = C
    for _ in range(20):
        x = [random.randrange(F.q) for _ in range(n)]
        vals = list(x)
        for a, u, b, v in gates:
            vals.append(F.add[F.mul[a][vals[u]]][F.mul[b][vals[v]]])
        for k, r in enumerate(rows):
            ref = 0
            for coef, xi in zip(r, x):
                ref = F.add[ref][F.mul[coef][xi]]
            if vals[outs[k]] != ref:
                return False
    return True


def improve_once(F, C, rows, DMAX, SMAX, WORKERS):
    n, gates, outs = C
    val = ws.values(F, C)
    probs = {}
    for S in ws.windows(C, SMAX):
        if not ws.is_convex(C, S):
            continue
        I, O = ws.window_problem(F, C, val, S)
        size = bin(S).count("1")
        if len(I) > DMAX or F.rank([val[i] for i in I]) < len(I):
            continue
        coords = [F.solve([val[i] for i in I], val[o]) for o in O]
        tg = sorted({F.normalize(tuple(c)) for c in coords})
        tg = [t for t in tg if sum(1 for x in t if x) > 1]
        wmax = size - 1 - len(tg)
        if wmax < 0:
            continue
        probs.setdefault((len(I), tuple(tg), wmax), (S, I))

    def job(item):
        (d, tg, wmax), (S, I) = item
        return item, ws.run_fs(d, F.q, list(tg), wmax, 20_000_000)

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for ((d, tg, wmax), (S, I)), (res, txt) in ex.map(job, probs.items()):
            if res != "FOUND":
                continue
            Wl = [l for l in txt.splitlines() if l.strip().startswith("W:")]
            Wp = [tuple(map(int, t.strip("()").split(","))) for t in Wl[0].split(":")[1].split()] if Wl else []

            def lift(c):
                v = tuple([0] * n)
                for ci, i in zip(c, I):
                    if ci:
                        v = F.lin(1, v, ci, val[i])
                return v
            outside = [val[n + g] for g in range(len(gates)) if not S >> (n + g) & 1]
            fin = {F.normalize(r): r for r in rows}
            newvals = [fin.get(F.normalize(v), v) for v in outside + [lift(c) for c in list(tg) + Wp]]
            C2 = rebuild(F, n, rows, newvals)
            if C2 and len(C2[1]) < len(gates) and check_random(F, C2, rows):
                ex.shutdown(wait=False, cancel_futures=True)
                return C2, len(probs)
    return None, len(probs)


def main():
    q, name = int(sys.argv[1]), sys.argv[2]
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    DMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    SMAX = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    WORKERS = int(sys.argv[6]) if len(sys.argv) > 6 else 12
    F = ws.GF(q)
    w = F.root(N)
    rows = [tuple(F.pw(w, k * m) for m in range(N)) for k in range(N)]
    C = base_circuits(F, N, w)[name]
    print(f"GF({q}) N={N}: start '{name}' with {len(C[1])} gates", flush=True)
    while True:
        C2, np_ = improve_once(F, C, rows, DMAX, SMAX, WORKERS)
        if C2 is None:
            print(f"  local minimum: {len(C[1])} gates ({np_} window problems, none improvable)", flush=True)
            break
        C = C2
        print(f"  -> {len(C[1])} gates (verified on random inputs)", flush=True)
    print("final gates (a, u, b, v):")
    for i, (a, u, b, v) in enumerate(C[1]):
        print(f"  n{N + i} = {a}*n{u} + {b}*n{v}")
    print("outputs X[k] = node", C[2])


if __name__ == "__main__":
    main()
