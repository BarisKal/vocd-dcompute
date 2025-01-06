# vocd-dcompute

Python implementation of computing D-optimum from VOCD.
For measuring lexical diversity (LD), a method called _D_ was developed by Malvern et al. The first version of the method was implemented by Brian MacWhinney in C++. Here, you can find a Python version of that implementation.

The original implementation can be found at https://dali.talkbank.org/clan/ (last visited 05.01.2025).

Detailed explanation of D can be found in:
Malvern D., Richards B., Chipere N., Durán P., 2004. Lexical diversity and language development: Quantification and assessment. Palgrave Macmillan UK.
https://link.springer.com/book/10.1057/9780230511804

---

# Example folder

In the example folder the output file "example*output_vocd.txt" was generated with the CHAT file "68.cha" (taken from /clan/examples/transcripts/ne32/68.cha).
The output was generated with CLANWin (v02.01.2025) with the command: vocd +t\*CHI +sm;*,o% 68.cha

D_optimum average in "example_output_vocd.txt" and the Python implementation give the same results.
