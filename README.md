# Parallel Longest Common Subsequence in Chapel

**What this is.** The longest common subsequence (LCS) of two strings is a classic dynamic-programming problem. It is used in diff tools, bioinformatics and plagiarism detection. This repository has a **Chapel** implementation of a parallel LCS algorithm that runs in **O(log³ n)** parallel time with **O(m·n)** processors (*m ≤ n* are the string lengths). It is the code behind our IEEE HPEC 2023 paper. The program also has a sequential baseline and per-phase timers (prefix sums, column minima, cell finding, crossing detection), so you can see where parallel time goes.

## Paper

S. Vahidi, B. Schieber, Z. Du, D. A. Bader. **Parallel Longest Common SubSequence Analysis in Chapel.** *2023 IEEE High Performance Extreme Computing Conference (HPEC).*
[doi:10.1109/HPEC58863.2023.10363472](https://doi.org/10.1109/HPEC58863.2023.10363472) · preprint [arXiv:2309.09072](https://arxiv.org/abs/2309.09072)

```bibtex
@inproceedings{vahidi2023parallel,
  author    = {Vahidi, Soroush and Schieber, Baruch and Du, Zhihui and Bader, David A.},
  title     = {Parallel Longest Common SubSequence Analysis In Chapel},
  booktitle = {2023 IEEE High Performance Extreme Computing Conference (HPEC)},
  year      = {2023},
  doi       = {10.1109/HPEC58863.2023.10363472}
}
```

## Files

- `lcs_project.chpl`: the whole program. For each of `numtest` test cases (600 by default) it reads the two strings from `newtest<k>.txt`, then times the parallel LCS (`paralcs`) and a sequential dynamic-programming LCS (`sequential`) on the same input. Average parallel and sequential times are printed per (log₂ n, log₂ m) size bucket. A results file (`168_8.txt` in the current code) is opened for output. Test-case files are not included in the repo.

## Build and run

Install [Chapel](https://chapel-lang.org/), then:

```bash
chpl --fast lcs_project.chpl -o lcs
./lcs
```

The program prints the number of available processing units (`here.numPUs()`, `here.maxTaskPar`) at startup. Put the `newtest<k>.txt` input files in the working directory first, with the first two lines holding the two strings. To change the number of tests or the output file, edit `numtest` and the `open(...)` call at the top of the file.

## License

MIT, see `LICENSE`.
