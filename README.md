# Parallel Longest Common Subsequence (Chapel)

A **Chapel** implementation of a **parallel algorithm** for the **longest common subsequence (LCS)** of two strings. The algorithm runs in **O(log^3 n)** time using **O(m n)** processors, where m and n are the lengths of the two strings (with n >= m).

## Reference

- **Paper (IEEE):** https://ieeexplore.ieee.org/abstract/document/10363472/

## What is in this repo

- Chapel source code that takes two strings (e.g. string1 and string2) and computes their LCS.
- Instructions to compile and run with the Chapel compiler.
- Example inputs or test scripts (if present).

## How to run

1. Install [Chapel](https://chapel-lang.org/) and set up the environment.
2. Compile the Chapel program (e.g. chpl lcs.chpl).
3. Run with two input strings (from file or command line, as documented in the repo).

## Complexity

- Time: O(log^3 n) parallel steps.
- Processors: O(m n).

## License

See the LICENSE file in the repository. For academic use, please cite the IEEE paper above and this repository.
