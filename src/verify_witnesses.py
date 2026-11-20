#!/usr/bin/env python3
"""Stand-alone check of the explicit circuits quoted in the SPL letter and in the doc.

No dependency on the search code, on Singular or on any package: exact arithmetic in
K = Z_p[y]/(m(y)) (p = 0 means Q, with Fractions).  For every circuit, each gate
z = a*u + b*v is expanded as a linear form in the inputs, and every claimed output is
compared with the corresponding row of the target matrix, exactly.

For the DFT circuits whose coefficients are all powers of omega, the script also counts
the edges of the circuit after the minimal conversion to an IFN in the sense of
Papadimitriou (J. ACM 1979): results must be sinks and inputs may not feed results
directly, so a relay node (one extra edge) is inserted after every output that feeds a
later gate, and after every input that feeds an output.

usage: python3 verify_witnesses.py
"""
from fractions import Fraction
from itertools import product


class Field:
    """K = Z_p[y]/(m), m monic, given by its coefficients m[0..deg-1] (y^deg = -sum m_j y^j)."""

    def __init__(self, p, m, name):
        self.p, self.m, self.d, self.name = p, m, len(m), name

    def red(self, c):
        return tuple(x % self.p for x in c) if self.p else tuple(Fraction(x) for x in c)

    def el(self, *c):
        c = list(c) + [0] * (self.d - len(c))
        return self.red(c)

    def add(self, a, b):
        return self.red([x + y for x, y in zip(a, b)])

    def mul(self, a, b):
        c = [0] * (2 * self.d - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i + j] += x * y
        for k in range(len(c) - 1, self.d - 1, -1):  # y^k = y^(k-d) * y^d
            for j in range(self.d):
                c[k - self.d + j] -= c[k] * self.m[j]
            c[k] = 0
        return self.red(c[:self.d])

    def pow(self, a, n):
        r = self.el(1)
        for _ in range(n):
            r = self.mul(r, a)
        return r

    def zero(self):
        return self.el(0)


def check(name, F, N, inputs, gates, outputs, rows, omega=None):
    """gates: list of (out_name, a, u, b, v); outputs: list of node names; rows: target rows."""
    form = {x: tuple(F.el(1) if j == i else F.zero() for j in range(N))
            for i, x in enumerate(inputs)}
    for z, a, u, b, v in gates:
        assert z not in form, z
        form[z] = tuple(F.add(F.mul(a, s), F.mul(b, t)) for s, t in zip(form[u], form[v]))
    ok = all(form[o] == tuple(r) for o, r in zip(outputs, rows))
    msg = f"{name}: {len(gates)} gates over {F.name}: outputs {'OK' if ok else 'WRONG'}"
    if omega is not None:
        powers = {F.pow(omega, k) for k in range(2 * N)}
        roots = all(a in powers and b in powers for _, a, _, b, _ in gates)
        outs = set(outputs)
        reused = {o for o in outs if any(o in (u, v) for _, _, u, _, v in gates)}
        direct = {x for x in inputs if any(z in outs and x in (u, v) for z, _, u, _, v in gates)}
        edges = 2 * len(gates) + len(reused) + len(direct)
        msg += (f"; multipliers are powers of omega: {roots}; IFN edges "
                f"{2 * len(gates)} + {len(reused)} reused outputs + {len(direct)} input relays "
                f"= {edges} (FFT: {2 * N * (N.bit_length() - 1)})")
    print(msg)
    return ok


def dft_rows(F, N, w):
    return [[F.pow(w, (k * n) % N) for n in range(N)] for k in range(N)]


def main():
    ok = True

    # 3-point DFT over Q(omega), omega^2 + omega + 1 = 0: 5 gates
    K = Field(0, [1, 1], "Q(omega), omega^2+omega+1=0")
    w, one = K.el(0, 1), K.el(1)
    m1 = K.el(-1)
    w2 = K.mul(w, w)
    ok &= check("F3", K, 3, ["x0", "x1", "x2"], [
        ("w0", one, "x1", one, "x2"),
        ("X0", one, "x0", one, "w0"),
        ("w1", one, "X0", K.add(w, m1), "x2"),
        ("X2", one, "w1", K.add(w2, m1), "x1"),
        ("X1", K.add(one, K.mul(m1, w)), "x0", w, "w1"),
    ], ["X0", "X1", "X2"], dft_rows(K, 3, w))

    # 5-point DFT over Q(omega), Phi_5: 15 gates, Lev (1980), Sec. 5.4.3 ("improved
    # algorithm"), with S2 = y1 - y4 and S4 = y3 - y2 as in Winograd's algorithm
    K5 = Field(0, [1, 1, 1, 1], "Q(omega), Phi_5(omega)=0")
    w = K5.el(0, 1)
    pw = [K5.pow(w, k) for k in range(5)]
    half = Fraction(1, 2)
    sc = lambda c, a: K5.mul(K5.el(c), a)
    g = sc(half, K5.add(pw[1], pw[4]))                       # cos u
    h = sc(half, K5.add(pw[2], pw[3]))                       # cos 2u
    f = sc(half, K5.add(pw[1], sc(-1, pw[4])))               # i sin u
    d = sc(half, K5.add(pw[2], sc(-1, pw[3])))               # i sin 2u
    hinv = sc(-2, K5.add(pw[1], pw[4]))                      # 1/h
    assert K5.mul(h, hinv) == K5.el(1)
    one, m1 = K5.el(1), K5.el(-1)
    ok &= check("F5 (Lev 1980)", K5, 5, ["y0", "y1", "y2", "y3", "y4"], [
        ("S1", one, "y1", one, "y4"),
        ("S2", one, "y1", m1, "y4"),
        ("S3", one, "y3", one, "y2"),
        ("S4", one, "y3", m1, "y2"),
        ("E1", one, "y0", h, "S1"),
        ("E2", one, "E1", h, "S3"),
        ("X0", hinv, "E2", K5.add(one, sc(-1, hinv)), "y0"),
        ("S10", one, "E2", K5.add(g, sc(-1, h)), "S1"),
        ("S11", one, "E1", g, "S3"),
        ("S12", f, "S2", sc(-1, d), "S4"),
        ("S13", d, "S2", f, "S4"),
        ("X1", one, "S10", one, "S12"),
        ("X4", one, "S10", m1, "S12"),
        ("X2", one, "S11", one, "S13"),
        ("X3", one, "S11", m1, "S13"),
    ], ["X0", "X1", "X2", "X3", "X4"], dft_rows(K5, 5, w))

    # 5-point DFT, Winograd's algorithm as written by Lev (Sec. 5.4.3): 17 additions,
    # multiplications M1..M5 absorbed into the gates
    a = K5.add(sc(half, K5.add(g, h)), m1)                   # (cos u + cos 2u)/2 - 1
    b = sc(half, K5.add(g, sc(-1, h)))                       # (cos u - cos 2u)/2
    c = K5.add(f, d)
    e = K5.add(f, sc(-1, d))
    ok &= check("F5 (Winograd 1978)", K5, 5, ["y0", "y1", "y2", "y3", "y4"], [
        ("S1", one, "y1", one, "y4"),
        ("S2", one, "y1", m1, "y4"),
        ("S3", one, "y3", one, "y2"),
        ("S4", one, "y3", m1, "y2"),
        ("S5", one, "S1", one, "S3"),
        ("S6", one, "S1", m1, "S3"),
        ("S7", one, "S2", one, "S4"),
        ("X0", one, "S5", one, "y0"),
        ("S9", one, "X0", a, "S5"),
        ("S10", one, "S9", b, "S6"),
        ("S11", one, "S9", sc(-1, b), "S6"),
        ("S12", c, "S2", sc(-1, d), "S7"),
        ("S13", d, "S7", e, "S4"),
        ("X1", one, "S10", one, "S12"),
        ("X4", one, "S10", m1, "S12"),
        ("X2", one, "S11", one, "S13"),
        ("X3", one, "S11", m1, "S13"),
    ], ["X0", "X1", "X2", "X3", "X4"], dft_rows(K5, 5, w))

    # 4-point Walsh-Hadamard over Q: 7 gates with a unit coefficient in each
    Q = Field(0, [0], "Q")
    one, m1, two, m2 = Q.el(1), Q.el(-1), Q.el(2), Q.el(-2)
    H4 = [[Q.el((-1) ** bin(k & n).count("1")) for n in range(4)] for k in range(4)]
    ok &= check("H4", Q, 4, ["x0", "x1", "x2", "x3"], [
        ("w0", one, "x2", one, "x3"),
        ("w1", one, "x1", one, "w0"),
        ("X0", one, "x0", one, "w1"),
        ("X2", one, "X0", m2, "w0"),
        ("w2", one, "X2", two, "x2"),
        ("X1", one, "w2", m2, "x1"),
        ("X3", two, "x0", m1, "w2"),
    ], ["X0", "X1", "X2", "X3"], H4)

    # 4-point DFT over F_5, omega = 2: 7 gates
    F5 = Field(5, [0], "F_5")
    one, m2 = F5.el(1), F5.el(-2)
    ok &= check("F4 (char 5)", F5, 4, ["x0", "x1", "x2", "x3"], [
        ("w0", one, "x1", one, "x3"),
        ("w1", one, "x2", one, "w0"),
        ("X0", one, "x0", one, "w1"),
        ("X2", one, "X0", m2, "w0"),
        ("w2", one, "X2", m2, "w1"),
        ("X1", one, "w2", one, "x3"),
        ("X3", one, "w2", one, "x1"),
    ], ["X0", "X1", "X2", "X3"], dft_rows(F5, 4, F5.el(2)), omega=F5.el(2))

    # 5-point DFT over F_11, omega = 3: 14 gates (local improvement of Lev's circuit;
    # wiring with integer obstruction d = 11, coefficients = the unique solution mod 11)
    F11 = Field(11, [0], "F_11")
    e = lambda c: F11.el(c)
    ok &= check("F5 (char 11)", F11, 5, ["x0", "x1", "x2", "x3", "x4"], [
        ("g1", e(1), "x1", e(1), "x4"),
        ("g2", e(1), "x1", e(-1), "x4"),
        ("g3", e(1), "x0", e(-4), "g1"),
        ("g4", e(1), "x2", e(1), "x3"),
        ("g5", e(1), "g2", e(4), "g4"),
        ("g6", e(1), "x3", e(-2), "g5"),
        ("g7", e(1), "g3", e(-4), "g4"),
        ("g8", e(-1), "g6", e(1), "g7"),
        ("g9", e(1), "g3", e(1), "g6"),
        ("g10", e(1), "x2", e(4), "g5"),
        ("g11", e(4), "x0", e(-3), "g7"),
        ("g12", e(1), "g1", e(-5), "g7"),
        ("g13", e(4), "g10", e(2), "g12"),
        ("g14", e(-4), "g10", e(2), "g12"),
    ], ["g11", "g13", "g8", "g9", "g14"], dft_rows(F11, 5, e(3)))

    # 8-point DFT over F_9 = F_3[i], i^2 = -1, omega = 1 + i: 23 gates
    F9 = Field(3, [1, 0], "F_9 = F_3[i]")
    c = {"1": F9.el(1), "-1": F9.el(-1), "i": F9.el(0, 1), "-i": F9.el(0, -1),
         "1-i": F9.el(1, -1), "1+i": F9.el(1, 1), "-1+i": F9.el(-1, 1), "-1-i": F9.el(-1, -1)}
    g = [("g1", "1", "x0", "1", "x2"), ("g2", "1", "g1", "1", "x4"),
         ("g3", "1", "x2", "-1", "x6"), ("g4", "1", "g2", "1", "x6"),
         ("g5", "1", "g2", "1", "g3"), ("g6", "1", "x1", "1", "x3"),
         ("g7", "1", "g6", "1", "x5"), ("g8", "1", "x3", "-1", "x7"),
         ("g9", "1", "g7", "1", "x7"), ("g10", "1", "g7", "1", "g8"),
         ("g11", "1-i", "g8", "1", "g9"), ("g12", "-1", "g6", "-1", "g11"),
         ("g13", "1", "g4", "1", "g9"), ("g14", "1", "g4", "-1", "g9"),
         ("g15", "1", "g5", "-i", "g10"), ("g16", "1", "g5", "i", "g10"),
         ("g17", "1", "g1", "1", "g4"), ("g18", "1", "g3", "-1+i", "g17"),
         ("g19", "1", "g12", "1", "g18"), ("g20", "1-i", "g8", "-1-i", "g19"),
         ("g21", "1+i", "g18", "-1", "g20"), ("g22", "-1-i", "g17", "-1+i", "g19"),
         ("g23", "-1+i", "g12", "1", "g22")]
    gates = [(z, c[a], u, c[b], v) for z, a, u, b, v in g]
    outs = ["g13", "g21", "g15", "g23", "g14", "g20", "g16", "g22"]
    ok &= check("F8 (char 3)", F9, 8, [f"x{j}" for j in range(8)], gates, outs,
                dft_rows(F9, 8, c["1+i"]), omega=c["1+i"])

    print("ALL OK" if ok else "FAILED")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
