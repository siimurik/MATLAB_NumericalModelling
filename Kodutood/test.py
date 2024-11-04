import sympy as sp

x = sp.symbols('x')
y = sp.Function('y')(x)
vorrand = x * y.diff(x) - y * (y - 1)
yldlahend = sp.dsolve(vorrand, y)

print(yldlahend)
