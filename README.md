# Exact additive complexity of small DFT kernels: code and logs

Code and computation logs for the letter *Exact Additive Complexity of Small DFT Kernels*
(A. Ugolini), its supplementary material and its information for reproducibility.

All the runs of the letter were made with the version of the code tagged `runs` (SHA-256 of
`groebner/exact_small.py`: `55718f8ad9157e5adbb28eb4e7a0b1c96909f721286453f65e59000ac673b896`).
Every log of `exact_small.py` contains a manifest line with the hash of the source, the
version of Singular and the options, and every answer of Singular is checked: an error
message or an unexpected output aborts the run instead of discarding a branch
(`tests/test_backend_failure.py`). The release `v1.0` adds checks and tests without changing
the enumeration. The Slurm outputs contain a line `code: <hash>` with the revision of the
private repository in which the runs were prepared; the SHA-256 of the source in the
manifest of every log identifies the code. The runs for F5 with G = 10, 11 and 12 were
executed twice, the first time with earlier versions of the code that did not
check the answers of Singular; the numbers of search nodes and of Gröbner bases coincide in
every part, and the logs of the first execution are in `results/first_execution/`. The environment is in
`environment.yml` (Python 3.11, Singular 4.4.1, SymPy, a C compiler).

The model counts the gates `z = a*u + b*v` of a linear circuit (one addition each),
with arbitrary constants `a, b` treated as free. The searches prove lower bounds
(no circuit with `G` gates exists) and find explicit circuits.

## Requirements

- Python 3.8 or later; `sympy` for `src/independent_check.py` and `groebner/obstruction.py`
- [Singular](https://www.singular.uni-kl.de) 4.4.1 for `groebner/exact_small.py` and
  `groebner/obstruction.py` (on `PATH`, or through the `SINGULAR` environment variable)
- a C compiler for `src/fs_search.c` (`make -C src`, GNU C99)

## Contents

| File | Purpose |
|---|---|
| `groebner/exact_small.py` | exhaustive search over circuit wirings with Gröbner basis tests over Q(ω); proves L > G over C |
| `groebner/obstruction.py` | integer obstruction d of a fixed wiring (Gröbner basis over Z) |
| `src/independent_check.py` | independent exhaustive check (sympy, no normalization, no ordering rules) of the small negative results |
| `src/check_certificates.py` | independent check of the negative results with certificates: own enumeration, no computer algebra system, exact verification of the identities sum h_i f_i = 1 written by `exact_small.py --cert` |
| `src/verify_witnesses.py` | exact check of every explicit circuit (F3, F5 of Lev and Winograd, H4, F4 over GF(5), F5 over GF(11), F8 over GF(9)) |
| `src/fs_search.c` | exact minimum over GF(p) or GF(p^2) (projective search) |
| `src/window_search.py`, `src/descent.py` | exact local improvement of circuits (23-gate F8 over GF(9); F5 from the circuits of Lev and Winograd) |
| `results/` | logs of all the runs reported in the letter |
| `tests/test_backend_failure.py` | regression test: a failure of Singular must abort the search, never become an exclusion |
| `tests/test_positive_control.py` | positive control: Lev's 15-gate circuit (and the 14-gate circuit over GF(11)), in the order used in the proof, survives the prime-length and Galois reductions and reaches a leaf |
| `cluster/rerun/` | Slurm scripts of the large runs (job-local memory watchdog) |
| `tools/summarize_parts.py` | checks a split run (every part complete, no error) and sums its statistics |
| `src/explore/cert_split.py`, `src/explore/check_certificates_split.py` | certificates for one part of a split run, and their independent check restricted to that part (wrappers: `exact_small.py` and `check_certificates.py` are imported, not modified) |
| `cluster/cert11/run_cert_part.sh` | Slurm script of the certificates for F5 with G = 11 (48 parts) |

## Results and logs

| Claim | Log |
|---|---|
| L(F3) = 5 over C (no circuit with 4 gates; 138 leaves with 5) | `results/exact_small.log`, `results/exact_small_depth.log` |
| L(F4) = 8 over C (no circuit with 7 gates) | `results/exact_dft4_G7.log` |
| L(H4) = 7 over C (no circuit with 6 gates; 1880 leaves with 7) | `results/exact_small.log` |
| L(F5) >= 11 over C (no circuit with 10 gates, 12 parts) | `results/exact_dft5_G10/` |
| L(F5) >= 12 over C (no circuit with 11 gates, 48 parts) | `results/exact_dft5_G11/` |
| L(F5) >= 13 over C (no circuit with 12 gates, 2048 parts) | `results/exact_dft5_G12/` (`summary.txt`) |
| L(F4) = 7 in characteristic 5 (no circuit with 6 gates; 16 leaves with 7) | `results/exact_dft4_G6_char5.log`, `results/exact_dft4_G7_char5.log` |
| Cross-checks: no ordering rules, no prime-length reduction, Galois reduction, independent implementation | `results/exact_small_depth.log`, `results/exact_dft5_G10_norules/`, `results/exact_dft5_G9_nosym2.log`, `results/exact_galois_small.log`, `results/exact_dft5_G10_galois/`, `results/independent_dft4_G7/` |
| Integer obstructions d = 5 (F4, char. 5), d = 3 (F8, char. 3), d = 11 (F5, char. 11) | `results/obstruction.log` |
| No circuit with 12 gates for F5 with coefficients in GF(11) (256 parts) | `results/fs_dft5_q11_w7/` |
| No window problem of Lev's 15-gate F5 circuit (connected convex windows of at most 10 gates, at most 5 distinct and linearly independent entering signals) admits one gate less over GF(31) and GF(41); a statement on these windows, not a lower bound (over GF(41), the job reached the 3-day time limit of the cluster after all the problems of Lev's circuit had ended, while some of Winograd's were still running: its log has no final summary) | `results/window_f5/` (`q31/`, `q41/`) |
| 23-gate F8 circuit over GF(9) | `results/descent_q9.log`, `results/win_q9.log` |
| All explicit circuits verified exactly | `results/verify_witnesses.log` |
| Certificates for every discarded branch: F3 (G = 4), H4 (G = 6), F4 (G = 7), F4 in characteristic 5 (G = 6), F5 (G = 9 and G = 10, i.e., L(F5) >= 11 without Singular), all checked | `results/certificates/` (`check.log`, `dft5_G10.log`) |
| Certificates for F5 with G = 11, i.e., L(F5) >= 12 without Singular: 2,199,741,770 certificates, generated and checked in 48 parts (about 1,110 core-hours); not stored (about 21 GB): the generation is deterministic, and the SHA-256 of each part, computed on the uncompressed stream, is published | `results/cert11/G11_logs/` (`gen*.log`, `check*.log`), `results/cert11_post/` (`sha*.txt`) |

## Reproducing the results

Run the searches from `groebner/`:

```sh
python3 exact_small.py dft3 4                    # no circuit: L(F3) > 4
python3 exact_small.py dft3 5                    # circuits exist (138 leaves)
python3 exact_small.py dft4 7                    # no circuit: L(F4) > 7
python3 exact_small.py wht4 6                    # no circuit: L(H4) > 6
python3 exact_small.py wht4 7                    # circuits exist (1880 leaves)
python3 exact_small.py dft4 6 --char 5 --omega 2 # no circuit in characteristic 5
python3 exact_small.py dft4 7 --char 5 --omega 2 # 16 leaves in characteristic 5
# F5: split into parts, k = 0..P-1, all parts must report "RESULT: no circuit"
python3 exact_small.py dft5 10 --sym2 --split k 12
python3 exact_small.py dft5 11 --sym2 --split k 48 --split-depth 4
# G = 12 with the Galois reduction (first target with a nonzero label fixed to row 1)
python3 exact_small.py dft5 12 --sym2 --galois --split k 2048 --split-depth 6
```

The Galois reduction (`--galois`, prime N, characteristic 0 only) is checked on
the smaller instances:

```sh
python3 exact_small.py dft3 5 --galois           # circuits exist (108 leaves instead of 138)
python3 exact_small.py dft5 10 --sym2 --galois --split k 12
```

Cross-checks (same conclusions, more Gröbner bases):

```sh
python3 exact_small.py dft3 4 --no-rules --no-support
python3 exact_small.py dft4 7 --no-rules
python3 exact_small.py dft5 9                    # without the prime-length reduction
python3 exact_small.py dft5 10 --sym2 --no-rules --split k 12
python3 ../src/independent_check.py dft3 4       # also dft3 5, wht4 6
python3 ../src/independent_check.py dft4 7 --split k 6
```

Certificates (the checker enumerates the sequences with its own code and discards a branch
only if its certificate verifies; about 1 minute for all the files except F5 with G = 10):

```sh
python3 exact_small.py dft4 7 --cert cert.jsonl.gz   # write the certificates (30 s)
python3 ../src/check_certificates.py dft4 7 ../results/certificates/dft4_G7.jsonl.gz
python3 ../src/check_certificates.py dft3 4 ../results/certificates/dft3_G4.jsonl.gz
python3 ../src/check_certificates.py wht4 6 ../results/certificates/wht4_G6.jsonl.gz
python3 ../src/check_certificates.py dft4 6 ../results/certificates/dft4_G6_char5.jsonl.gz --char 5 --omega 2
python3 ../src/check_certificates.py dft5 9 ../results/certificates/dft5_G9_sym2.jsonl.gz --sym2
# F5 with G = 10: 12,780,390 certificates in 3 files (104 MB), read in order (5.4 h, 14 MB)
python3 ../src/check_certificates.py dft5 10 '../results/certificates/dft5_G10_sym2.part*.jsonl.gz' --sym2 --stream
# F5 with G = 11: part k of 48 (k = 0..47), with exact_small.py and check_certificates.py of tag
# `runs` (git checkout runs -- groebner/exact_small.py src/check_certificates.py), about 23
# core-hours per part; the hash must equal results/cert11_post/sha$k.txt
python3 ../src/explore/cert_split.py dft5 11 cert$k.jsonl.gz --split $k 48 --split-depth 4 --sym2
python3 ../src/explore/check_certificates_split.py dft5 11 cert$k.jsonl.gz --split $k 48 --split-depth 4 --sym2 --stream
gzip -dc cert$k.jsonl.gz | sha256sum
```

Explicit circuits and integer obstructions:

```sh
python3 ../src/verify_witnesses.py               # prints ALL OK
python3 ../tests/test_backend_failure.py         # prints PASS (failure handling)
python3 ../tests/test_positive_control.py        # prints PASS (known circuits survive the reductions)
python3 obstruction.py f4_char5                  # d = 5
python3 obstruction.py f8_char3                  # d = 3
python3 obstruction.py f5_char11                 # d = 11 (14-gate F5 over GF(11))
```

Searches over finite fields, run from `src/` (`make` first). The program reports
L = r + |W|, with W the auxiliary points of the circuit; two quick checks:

```sh
./fs_search 4 5 0 3 dft                          # F4 over GF(5): none with L <= 6, found with L = 7
./fs_search 3 7 0 2 dft                          # F3 over GF(7): found with L = 5
# exact search over GF(11): no circuit with 12 gates (256 parts, k = 0..255)
./fs_search 5 11 7 7 dft split k 256 all 1000   # all N: print up to N circuits
# local improvement of the circuits of Lev (15 gates) and Winograd (17) over GF(q)
python3 window_search.py 31 5 5 10 2e9 16        # also q = 41; node limit 2e9 per problem, 16 threads
DUMP=dir python3 window_search.py 31 5 5 10      # write the window problems only
./fs_search 5 31 4 4 file:dir/q31_probNNNN.txt maxnodes 2000000000
```

The searches for F3, F4 and H4 with `exact_small.py`, the explicit circuits, the integer
obstructions and the two quick checks of `fs_search` run in less than a minute each on a
laptop, with less than 100 MB of memory; the independent check of F4 takes about 45 minutes
(6 parts), and the runs of F5 take from hours to thousands of CPU hours (see the
information for reproducibility of the letter). Every search log contains the
manifest, the statistics line printed by the program and a final `RESULT:` line; a part
without a `RESULT:` line did not complete. The runs of F5 with G >= 10 were split over
processes by a counter of the moves at a fixed depth (`--split k P --split-depth d`), and `tools/summarize_parts.py` checks every
part, the manifests and that the parts form a partition, and sums the CPU times (user and
system, including Singular). Any single part can be repeated and compared with its log.

## Citation and license

Citation: to be added after publication.

License: the code is released under the MIT License (`LICENSE`); the logs in
`results/` are released under the Creative Commons Attribution 4.0 International
License (CC BY 4.0, `results/LICENSE`).
