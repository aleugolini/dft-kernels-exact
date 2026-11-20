#!/usr/bin/env python3
"""Step 1: in which characteristics can a given circuit *wiring* compute the DFT?

A wiring fixes, for every gate g, the two nodes it combines (u_g, v_g) and which
gates are the outputs X[k].  The coefficients become unknowns a_g, b_g, omega is a
root of the N-th cyclotomic polynomial, and "the circuit computes F_N" is a system
of polynomial equations with integer coefficients.  A Groebner basis over Z gives
the ideal I; its intersection with Z is (d):
  * d = 0 and I != (1) over Q  -> the wiring works over C (characteristic 0);
  * d != 0                      -> it can only work in characteristics p | d.
Gauge fixing: every gate value may be rescaled (its consumers absorb the factor),
so a_g = 1 for gates that are not outputs (this assumes a_g != 0, i.e. that no
gate is a mere copy of a scaled node, which would make the circuit shorter).

usage: python3 obstruction.py CIRCUIT [--no-gauge] [--keep]
       CIRCUIT in: f4_char5, f8_char3, f5_char11  (see CIRCUITS below)
Requires Singular (brew install singular) and sympy.
"""
import os
import subprocess
import sys
import tempfile

import sympy as sp

# node numbering: 0..N-1 inputs, N.. gates in order; outputs[k] = node of X[k]
CIRCUITS = {
    # 7-gate 4-point DFT found over GF(5) (omega = 2)
    "f4_char5": dict(N=4, gates=[(1, 3), (2, 4), (0, 5), (6, 4), (5, 7), (3, 8), (1, 8)],
                     outputs=[6, 9, 7, 10]),
    # 14-gate 5-point DFT over GF(11) (omega = 3), local improvement of Lev's 15-gate circuit
    "f5_char11": dict(N=5, gates=[(1, 4), (1, 4), (0, 5), (2, 3), (6, 8), (3, 9), (7, 8), (10, 11),
                                  (7, 10), (2, 9), (0, 11), (5, 11), (14, 16), (14, 16)],
                      outputs=[15, 17, 12, 13, 18]),
    # 23-gate 8-point DFT found over GF(9) (omega = 1+i), results/descent_q9.log
    "f8_char3": dict(N=8, gates=[(0, 2), (4, 8), (2, 6), (6, 9), (9, 10), (1, 3), (5, 13), (3, 7),
                                 (7, 14), (14, 15), (15, 16), (13, 18), (11, 16), (11, 16),
                                 (12, 17), (12, 17), (8, 11), (10, 24), (19, 25), (15, 26),
                                 (25, 27), (24, 26), (19, 29)],
                     outputs=[20, 28, 22, 30, 21, 27, 23, 29]),
}


def system(c, gauge=True):
    N, gates, outs = c["N"], c["gates"], c["outputs"]
    w = sp.Symbol("w")
    A = sp.symbols(f"a1:{len(gates) + 1}")
    B = sp.symbols(f"b1:{len(gates) + 1}")
    out_nodes = set(outs)
    val = [[sp.Integer(1 if j == i else 0) for j in range(N)] for i in range(N)]
    subs = {}
    for g, (u, v) in enumerate(gates):
        node = N + g
        a = A[g]
        if gauge and node not in out_nodes:
            subs[a] = 1
            a = sp.Integer(1)
        val.append([sp.expand(a * val[u][j] + B[g] * val[v][j]) for j in range(N)])
    eqs = []
    for k, node in enumerate(outs):
        for n in range(N):
            target = w ** ((k * n) % N)
            e = sp.expand(val[node][n] - target)
            if e != 0:
                eqs.append(e)
    phi = sp.cyclotomic_poly(N, w)
    eqs.append(sp.expand(phi))
    vars_ = [w] + [a for a in A if a not in subs] + list(B)
    return eqs, vars_


def singular(eqs, vars_, coeff):
    names = ",".join(str(v) for v in vars_)
    polys = ",\n  ".join(str(e).replace("**", "^") for e in eqs)
    script = f"""option(redSB);
ring R = ({coeff}), ({names}), dp;
ideal I =
  {polys};
ideal G = std(I);
int i; int found = 0;
for (i = 1; i <= size(G); i++) {{ if (deg(G[i]) == 0) {{ print("CONSTANT " + string(G[i])); found = 1; }} }}
print("SIZE " + string(size(G)));
print("DIM " + string(dim(G)));
quit;
"""
    with tempfile.NamedTemporaryFile("w", suffix=".sing", delete=False) as f:
        f.write(script)
        path = f.name
    out = subprocess.run(["Singular", "-q", path], capture_output=True, text=True).stdout
    os.unlink(path)
    return out


def main():
    name = sys.argv[1]
    gauge = "--no-gauge" not in sys.argv
    eqs, vars_ = system(CIRCUITS[name], gauge)
    print(f"{name}: {len(eqs)} equations, {len(vars_)} unknowns, gauge={'on' if gauge else 'off'}")
    print("degrees:", sorted({sp.Poly(e, *vars_).total_degree() for e in eqs}))
    for coeff in ["integer", "0"]:
        out = singular(eqs, vars_, coeff)
        label = "over Z" if coeff == "integer" else "over Q"
        print(f"--- {label}:")
        print("   " + "\n   ".join(l for l in out.splitlines() if l.strip()))


if __name__ == "__main__":
    main()
