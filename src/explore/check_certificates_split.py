#!/usr/bin/env python3
"""Independent check of ONE PART of a split run with certificates (cert_split.py).

src/check_certificates.py (the letter's checker: own enumeration, no computer algebra system) is
read and executed with two filters inserted, which reproduce the split of exact_small.py: the moves
generated at gate depth SPLIT_DEPTH (targets first, then auxiliary nodes, in the order of the
search) are numbered in the order of the enumeration, and only those with number = k mod K are
followed; all the other moves are followed as in the unsplit check. The file check_certificates.py
is not modified. Running the K parts checks every branch of the tree exactly once; the per-part
counts must equal those of the search (Groebner tests = target moves of the checker).

usage: python3 check_certificates_split.py KERNEL G CERTFILE --split k K --split-depth d [--sym2] [--stream]
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CC = os.path.join(HERE, "..", "check_certificates.py")
src = open(CC).read()
T_HOOK = '                    stats["target_moves"] += 1\n'
W_HOOK = "                    form = [ring.add(nodes[u][3][c], ring.mul(t, nodes[v][3][c])) for c in range(N)]\n"
FILTER = ("                    if nG + 1 == SPLIT[2]:  # split (check_certificates_split.py)\n"
          "                        SPLIT[3][0] += 1\n"
          "                        if SPLIT[3][0] % SPLIT[1] != SPLIT[0]:\n"
          "                            continue\n")
if src.count(T_HOOK) != 1 or src.count(W_HOOK) != 1:
    raise SystemExit("hook lines not found in check_certificates.py")
src = src.replace(T_HOOK, FILTER + T_HOOK).replace(W_HOOK, FILTER + W_HOOK)
i = sys.argv.index("--split")
k, K = int(sys.argv[i + 1]), int(sys.argv[i + 2])
d = int(sys.argv[sys.argv.index("--split-depth") + 1])
ns = {"__name__": "check_certificates_split", "SPLIT": (k, K, d, [-1])}
exec(compile(src, CC, "exec"), ns)
print(f"part {k} of {K}, split depth {d}")
ns["main"]()
