"""
Символьный расчет одночастичной статистической суммы.
"""
import sympy as sp

p = sp.Symbol("p")
px, py, pz = sp.symbols("p_x p_y p_z")
m = sp.Symbol("\\left<m\\right>")

Z1 = 4*sp.pi*sp.sinh(p)/p

subs_exp_to_m = {
    sp.exp(2*p): (m*p+p+1)/(m*p-p+1),
    sp.cosh(p)/sp.sinh(p): m+1/p
}