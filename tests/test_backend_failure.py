#!/usr/bin/env python3
"""Regression test: a failure of Singular must abort the search, never become an exclusion.

On the known-positive instance F3 with G = 5 (138 leaves), the Groebner test of every target
move is corrupted with an undefined symbol, so Singular stays alive and reports an error.
The only acceptable outcome is an exception (aborted run); returning a result, in
particular 'no circuit', is a failure of the test. The same is checked for a corrupted
state command (auxiliary-node definition) and for the search over Z.

usage: python3 tests/test_backend_failure.py      (prints PASS or FAIL, exit status 0 on PASS)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "groebner"))
import exact_small as es  # noqa: E402


def corrupt(method, needle):
    orig = getattr(es.Singular, method)

    def bad(self, cmd):
        if needle in cmd:
            cmd = cmd.replace(needle, needle + "+undefinedsym_", 1)
        return orig(self, cmd)
    return orig, bad


def case(name, method, needle, **kw):
    orig, bad = corrupt(method, needle)
    setattr(es.Singular, method, bad)
    try:
        stats, _ = es.search("dft3", 5, show=0, **kw)
        print(f"FAIL {name}: the search returned {stats} instead of aborting")
        return False
    except RuntimeError as e:
        print(f"PASS {name}: aborted ({str(e).splitlines()[0]})")
        return True
    finally:
        setattr(es.Singular, method, orig)


def main():
    os.environ["PROGRESS_EVERY"] = "0"
    stats, _ = es.search("dft3", 5, show=0)
    ok = stats["leaves"] == 138
    print(f"{'PASS' if ok else 'FAIL'} control: F3, G=5 has {stats['leaves']} leaves (expected 138)")
    ok &= case("corrupted Groebner test", "run", "ideal(a")
    ok &= case("corrupted auxiliary node", "push", "+t1*")
    ok &= case("corrupted Groebner test over Z", "run", "ideal(a", overz=True)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
