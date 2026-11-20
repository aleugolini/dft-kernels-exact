#!/usr/bin/env python3
"""Positive control for the reductions used with N = 5: a known circuit must survive them.

For Lev's 15-gate circuit for F5 over C, the program builds, as in the proof of completeness,
1. the affine relabeling of the inputs that makes the first node x0 + t x1 (prime-length rule),
2. the greedy order (smallest key first, ties broken by a fixed numbering of the gates),
3. the Galois conjugation that makes the first target with a nonzero label the target 1,
and then runs exact_small.py with G = 15, --sym2 and --galois, restricted to this single
sequence (option only_path): every rule, the support test and the Groebner tests are active, and the
sequence must end in a leaf. The same is done in characteristic 11 for the 14-gate circuit over
F_11 (omega = 3), without the Galois rule, which is not used in positive characteristic.

usage: python3 tests/test_positive_control.py   (prints PASS or FAIL, exit status 0 on PASS)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "groebner"))
import exact_small as es  # noqa: E402

N = 5
LEV = [  # (gate, parent, parent); only the wiring matters
    ("S1", "y1", "y4"), ("S2", "y1", "y4"), ("S3", "y3", "y2"), ("S4", "y3", "y2"),
    ("E1", "y0", "S1"), ("E2", "E1", "S3"), ("X0", "E2", "y0"), ("S10", "E2", "S1"),
    ("S11", "E1", "S3"), ("S12", "S2", "S4"), ("S13", "S2", "S4"), ("X1", "S10", "S12"),
    ("X4", "S10", "S12"), ("X2", "S11", "S13"), ("X3", "S11", "S13"),
]
LEV_OUT = {"X0": 0, "X1": 1, "X2": 2, "X3": 3, "X4": 4}  # gate -> row k
CHAR11 = [
    ("g1", "y1", "y4"), ("g2", "y1", "y4"), ("g3", "y0", "g1"), ("g4", "y2", "y3"),
    ("g5", "g2", "g4"), ("g6", "y3", "g5"), ("g7", "g3", "g4"), ("g8", "g6", "g7"),
    ("g9", "g3", "g6"), ("g10", "y2", "g5"), ("g11", "y0", "g7"), ("g12", "g1", "g7"),
    ("g13", "g10", "g12"), ("g14", "g10", "g12"),
]
CHAR11_OUT = {"g11": 0, "g13": 1, "g8": 2, "g9": 3, "g14": 4}


def greedy(gates, out):
    """Greedy order of the nodes: at each position the available node with the smallest key,
    targets (T, k) before auxiliary nodes (W, (u, v)); ties broken by the numbering (list order)."""
    pos = {f"y{n}": n for n in range(N)}
    seq, left = [], list(range(len(gates)))
    while left:
        cand = []
        for i in left:
            g, a, b = gates[i]
            if a in pos and b in pos:
                u, v = sorted((pos[a], pos[b]))
                key = ("T", out[g]) if g in out else ("W", (u, v))
                cand.append((key, i, u, v))
        key, i, u, v = min(cand)
        g = gates[i][0]
        seq.append(("T", out[g], u, v) if g in out else ("W", None, u, v))
        pos[g] = N + len(seq) - 1
        left.remove(i)
    return seq


def relabel_inputs(gates, out):
    """Affine map n -> a n + b (mod N) sending the input pair of the first greedy node to (0, 1);
    row k becomes row k a^{-1} (up to a scale factor)."""
    first = greedy(gates, out)[0]
    assert first[0] == "W", "first node must be auxiliary for prime N"
    u, v = first[2], first[3]
    a = pow(v - u, -1, N)
    b = (-a * u) % N
    ren = lambda x: f"y{(a * int(x[1:]) + b) % N}" if x.startswith("y") else x
    g2 = [(g, ren(p), ren(q)) for g, p, q in gates]
    ainv = pow(a, -1, N)
    return g2, {g: (k * ainv) % N for g, k in out.items()}


def conjugate(gates, out):
    """Galois conjugation j = k0^{-1}, with k0 the smallest nonzero label available at the first
    target with a nonzero label, so that this target becomes the target 1."""
    seq = greedy(gates, out)
    pos = {f"y{n}": n for n in range(N)}
    for idx, node in enumerate(seq):
        if node[0] == "T" and node[1] != 0:
            k0 = node[1]
            break
    j = pow(k0, -1, N)
    return {g: (k * j) % N for g, k in out.items()}


def run(name, gates, out, G, **kw):
    gates, out = relabel_inputs(gates, out)
    if kw.get("galois"):
        out = conjugate(gates, out)
    path = greedy(gates, out)
    first_t = next(m for m in path if m[0] == "T" and m[1] != 0)
    stats, _ = es.search("dft5", G, show=0, sym2=True, only_path=path, **kw)
    ok = stats["leaves"] == 1 and path[0] == ("W", None, 0, 1) and (not kw.get("galois") or first_t[1] == 1)
    print(f"{'PASS' if ok else 'FAIL'} {name}: path of {len(path)} nodes, first node {path[0]}, "
          f"first nonzero target {first_t[1]}, leaves {stats['leaves']}, Groebner tests {stats['gb']}")
    return ok


def main():
    os.environ["PROGRESS_EVERY"] = "0"
    ok = run("Lev's 15-gate circuit over Q(omega), --sym2 --galois", LEV, LEV_OUT, 15, galois=True)
    ok &= run("Lev's 15-gate circuit over Q(omega), --sym2", LEV, LEV_OUT, 15)
    ok &= run("14-gate circuit over F_11 (omega = 3), --sym2", CHAR11, CHAR11_OUT, 14, char=11, omega=3)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
