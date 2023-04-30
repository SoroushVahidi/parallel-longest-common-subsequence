# parallel-longest-common-subsequence
It is a code in Chapel. Its input is 2 strings named string1 and string2, and it writes the lcs of them in O(log^3(n)) with mn processors, such that m and n are the lengths of the input strings and n>=m. You should set the input in the code. None of string1 and string2 can be empty. The algorithm is based on https://ieeexplore.ieee.org/document/298210
