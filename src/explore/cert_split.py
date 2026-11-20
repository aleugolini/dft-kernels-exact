#!/usr/bin/env python3
"""One part of a split run of groebner/exact_small.py WITH certificates.

exact_small.py refuses --cert together with --split in its command line, but its function search()
supports both: this wrapper calls search() directly (the file is imported, not modified), writes the
certificates of the part to CERTFILE (gzip if it ends in .gz), and prints the same lines as
exact_small.py (manifest, statistics, per-depth counts, RESULT), so that the run can be compared
part by part with the run without certificates (same split, same counts of Groebner bases).

usage: python3 cert_split.py KERNEL G CERTFILE --split k K --split-depth d [--sym2]
"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ES_PATH = os.path.join(HERE, "..", "..", "groebner", "exact_small.py")
sys.path.insert(0, os.path.dirname(ES_PATH))
import exact_small as es  # noqa: E402

kernel, G, cert = sys.argv[1], int(sys.argv[2]), sys.argv[3]
i = sys.argv.index("--split")
split = (int(sys.argv[i + 1]), int(sys.argv[i + 2]))
sd = int(sys.argv[sys.argv.index("--split-depth") + 1])
sym2 = "--sym2" in sys.argv
manifest = dict(source=os.path.basename(ES_PATH), source_sha256=hashlib.sha256(open(ES_PATH, "rb").read()).hexdigest(),
                argv=sys.argv[1:], split_depth=sd, python=sys.version.split()[0], wrapper="cert_split.py")
stats, sols = es.search(kernel, G, 0, None, 0, sym2, split, sd, True, True, False, False, cert, manifest)
print("manifest:", json.dumps(manifest, sort_keys=True))
per_depth = stats.pop("depth")
stats.pop("obstr")
print(f"{kernel}, at most {G} gates, over Q(omega)-bar, sym2={sym2}, split={split}: {stats}")
print("per depth (gate: moves, Groebner tests, rejected):",
      "; ".join(f"{d}: {m}, {g}, {r}" for d, (m, g, r) in per_depth.items()))
print("RESULT:", "circuits exist" if stats["leaves"] else f"no circuit -> L > {G}")
