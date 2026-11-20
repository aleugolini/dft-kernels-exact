#!/usr/bin/env python3
"""Summarize a split run of groebner/exact_small.py (one log per part) for the papers.

Checks that every part 0..P-1 has a log ending with 'RESULT: no circuit', with the expected
header (kernel, gates, split), no traceback/error, and no WATCHDOG line in the Slurm outputs;
then sums the statistics and prints them in LaTeX form (thin-space thousands).

usage: python3 summarize_parts.py DIR P [--gates G] [--kernel dft5]
       DIR holds part{k}.log (and optionally slurm_*.out); P is the number of parts.
Exit status 0 only if the run is complete and clean.
Also importable: summarize(DIR, P, G, kernel) returns the totals as a dict.
"""
import ast
import glob
import json
import os
import re
import sys


def tex(n):
    s = f"{n:,}".replace(",", r"\,")
    return f"${s}$"


def summarize(d, P, G=None, kernel="dft5"):
    tot = dict(dfs=0, gb=0, leaves=0, sec=0.0, cpu=0.0, restarts=0)
    configs = {}  # run configuration from the manifest line (all but the part index) -> parts
    elapsed, rss, per_part_sec = [], [], []
    depth = {}
    issues = dict(missing=[], incomplete=[], bad=[], found=[], wrong_header=[], watchdog=[], mixed=[])
    for k in range(P):
        p = os.path.join(d, f"part{k}.log")
        if not os.path.exists(p):
            issues["missing"].append(k)
            continue
        txt = open(p, errors="replace").read()
        if re.search(r"Traceback|Error|error:|Singular died", txt):
            issues["bad"].append(k)
        res = re.findall(r"^RESULT: (.*)$", txt, re.M)
        if not res:
            issues["incomplete"].append(k)
            continue
        if not res[-1].startswith("no circuit"):
            issues["found"].append((k, res[-1]))
        m = re.search(r"^(\w+), at most (\d+) gates, .*split=\((\d+), (\d+)\).*?: (\{.*\})$", txt, re.M)
        if not m or m.group(1) != kernel or (G and int(m.group(2)) != G) \
                or int(m.group(3)) != k or int(m.group(4)) != P:
            issues["wrong_header"].append(k)
            continue
        st = ast.literal_eval(m.group(5))
        mm = re.search(r"^manifest: (\{.*\})$", txt, re.M)
        if mm:
            man = json.loads(mm.group(1))
            argv = list(man.get("argv", []))
            if "--split" in argv:
                argv[argv.index("--split") + 1] = "k"
            key = json.dumps(dict(man, argv=argv), sort_keys=True)
            configs.setdefault(key, []).append(k)
        for key in tot:
            tot[key] += st.get(key, 0)
        per_part_sec.append(st.get("sec", 0.0))
        dm = re.search(r"^per depth \(.*?\): (.*)$", txt, re.M)
        if dm:
            for item in dm.group(1).split(";"):
                g, vals = item.split(":")
                a = [int(x) for x in vals.split(",")]
                acc = depth.setdefault(int(g), [0, 0, 0])
                for i in range(3):
                    acc[i] += a[i]
        tm = re.search(r"maxrss=(\d+)KB elapsed=([\d.]+)s", txt)
        if tm:
            rss.append(int(tm.group(1)))
            elapsed.append(float(tm.group(2)))
    models = {}
    for s in sorted(glob.glob(os.path.join(d, "slurm_*.out"))):
        t = open(s, errors="replace").read()
        if "WATCHDOG" in t:
            issues["watchdog"].append(os.path.basename(s))
        m = re.search(r"Model name:\s*(.*)", t)
        if m:
            name = re.sub(r"\(R\)|\(TM\)|CPU|@.*|\s+Processor", "", m.group(1))
            name = re.sub(r"\s+", " ", name).strip()
            models[name] = models.get(name, 0) + 1
    if len(configs) > 1:  # parts run with different code, options or Singular versions
        issues["mixed"] = [sorted(v)[:5] for v in configs.values()]
    n = P - len(issues["missing"]) - len(issues["incomplete"]) - len(issues["wrong_header"])
    clean = not any(issues.values()) and n == P
    return dict(dir=d, P=P, n=n, clean=clean, issues=issues, tot=tot, hours=tot["sec"] / 3600,
                per_part_sec=per_part_sec, elapsed=elapsed, rss=rss, models=models, depth=depth,
                configs=[json.loads(c) for c in configs])


def report(r):
    tot, h = r["tot"], r["hours"]
    lines = [f"directory: {r['dir']}",
             f"parts complete: {r['n']}/{r['P']}" + ("" if r["clean"] else "   ** NOT CLEAN **")]
    labels = dict(missing="missing", incomplete="no RESULT line (running or failed)",
                  bad="traceback/error", found="RESULT other than no circuit",
                  wrong_header="unexpected header", watchdog="WATCHDOG in",
                  mixed="different run configurations (first parts of each)")
    for key, lst in r["issues"].items():
        if lst:
            lines.append(f"  {labels[key]}: {len(lst)}: {lst[:20]}{' ...' if len(lst) > 20 else ''}")
    lines += [f"search nodes (dfs):  {tot['dfs']:>18,}   {tex(tot['dfs'])}",
              f"Groebner bases:      {tot['gb']:>18,}   {tex(tot['gb'])}",
              f"leaves:              {tot['leaves']:>18,}",
              f"Singular restarts:   {tot['restarts']:>18,}   {tex(tot['restarts'])}",
              f"search time:         {h:>18,.1f} h (sum of 'sec' over parts, elapsed)"]
    if tot["cpu"]:
        lines.append(f"CPU time:            {tot['cpu'] / 3600:>18,.1f} h (user + system, with Singular)")
    for c in r.get("configs", []):
        lines.append("configuration: " + json.dumps(c, sort_keys=True))
    pp = r["per_part_sec"]
    if pp:
        lines.append(f"  per part:          {min(pp) / 3600:.2f} to {max(pp) / 3600:.2f} h, "
                     f"mean {sum(pp) / len(pp) / 3600:.2f} h")
    if r["elapsed"]:
        lines.append(f"elapsed (/usr/bin/time): {sum(r['elapsed']) / 3600:,.1f} h summed over "
                     f"{len(r['elapsed'])} parts; max resident memory {max(r['rss']) / 1024:.0f} MB")
    if r["models"]:
        lines.append("CPU models (Slurm jobs): " + "; ".join(f"{m} ({c})" for m, c in sorted(r["models"].items())))
    if r["depth"]:
        lines.append("per depth (gate: moves, Groebner tests, rejected), summed over parts:")
        lines.append("  " + "; ".join(f"{g}: {a[0]}, {a[1]}, {a[2]}" for g, a in sorted(r["depth"].items())))
    lines += ["", "LaTeX placeholders:",
              f"  Table I, Groebner bases:  {tex(tot['gb'])}",
              f"  CPU hours (rounded):      ${round(h, -2):,.0f}$".replace(",", r"\,"),
              f"  run record: {tex(tot['dfs'])} search nodes, {tex(tot['gb'])} Groebner bases, "
              f"{tex(tot['restarts'])} restarts, {h:,.0f} hours"]
    return "\n".join(lines)


def main():
    d, P = sys.argv[1], int(sys.argv[2])
    G = int(sys.argv[sys.argv.index("--gates") + 1]) if "--gates" in sys.argv else None
    kernel = sys.argv[sys.argv.index("--kernel") + 1] if "--kernel" in sys.argv else "dft5"
    r = summarize(d, P, G, kernel)
    print(report(r))
    sys.exit(0 if r["clean"] else 1)


if __name__ == "__main__":
    main()
