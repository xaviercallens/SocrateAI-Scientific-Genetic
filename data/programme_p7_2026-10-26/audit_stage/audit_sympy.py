"""P7 auditor: sympy checks of the non-trivial algebra in the human proofs (Tier X support for a Tier C review)."""
import sympy as sp

r, i, N, j, n = sp.symbols("r i N j n", positive=True)
p = r * i / (r * i + N - i) * (N - i) / N
q = (N - i) / (r * i + N - i) * i / N
print("p - r q =", sp.simplify(p - r * q))
print("stay identity residual =", sp.simplify(1 - p - q - (r * i**2 + (N - i) ** 2) / ((r * i + N - i) * N)))
# KU form equality (PG-C-0004)
pku = r / (1 + (r - 1) * j / n) * (j / n) * (n - j) / n
qku = 1 / (1 + (r - 1) * j / n) * ((n - j) / n) * (j / n)
print("KU p+ vs moranUp residual =", sp.simplify(pku - p.subs({i: j, N: n})))
print("KU p- vs moranDown residual =", sp.simplify(qku - q.subs({i: j, N: n})))
# first-step equation for h(i) = S_i/S_N, S_m = (1 - r^-m)/(1 - r^-1)
S = lambda m: (1 - r**(-m)) / (1 - 1 / r)
h = lambda m: S(m) / S(N)
print("first-step residual =", sp.simplify(p * h(i + 1) + q * h(i - 1) + (1 - p - q) * h(i) - h(i)))
# Price five-sum identity
M, A, B, C, D = sp.symbols("M A B C D", nonzero=True)
lhs = A / M * (D / A - B / M)
cov = (C - A / M * B - B / M * A + A / M * (B / M) * M) / M
rhs = cov + (D - C) / M
print("Price scalar residual =", sp.simplify(lhs - rhs))
print("Cov raw form residual =", sp.simplify(cov - (C - A * B / M) / M))
# Neutral limit of closed form
x = sp.symbols("x", positive=True)
for ii, NN in [(1, 3), (2, 5), (7, 10)]:
    print("limit r->1", ii, NN, sp.limit((1 - x**(-ii)) / (1 - x**(-NN)), x, 1), "expected", sp.Rational(ii, NN))
