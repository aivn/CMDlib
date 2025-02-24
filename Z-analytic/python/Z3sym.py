"""
Символьный расчет трехчастичной статистической суммы.
"""
import sympy as sp
import symbase
import Z2sym
import Z3num

SIGMAS = [-1, 1]

hi, hj, lij, lik = sp.symbols("h_i h_j \\lambda_{ij} \\lambda_{ik}", positive=True)
s1, s2, s3 = symbase.SigmaSymbol("\\sigma_1"), symbase.SigmaSymbol("\\sigma_2"), symbase.SigmaSymbol("\\sigma_3")

mi = sp.Symbol("\\left<{m_i}\\right>")
mj = sp.Symbol("\\left<{m_j}\\right>")
eta_ij = sp.Symbol("\\left<{\\eta_{ij}}\\right>")
eta_ik = sp.Symbol("\\left<{\\eta_{ik}}\\right>")

eta_ik_2 = sp.Symbol("\\left<\eta_{ik}^{2}\\right>")
mj_eta_ij = sp.Symbol("\\left<m_j\\eta_{ij}\\right>")
mi_eta_ik = sp.Symbol("\\left<m_i\\eta_{ik}\\right>")
mj_eta_ik = sp.Symbol("\\left<m_j\\eta_{ik}\\right>")
mjh2 = sp.Symbol("\\left<m_{jh}^2\\right>")
Q_star = sp.Symbol("Q^*")


# ===================================== Модули для sp.lambdify =====================================
num_module_norm = Z2sym.num_module_norm[0]
num_module_norm.update({
    "EFunction": lambda *args: Z3num.E_norm_symmetrical(*args[:-3], *map(bool, args[-3:])),
    "FFunction": lambda *args: Z3num.F_norm_symmetrical(*args[:-3], *map(bool, args[-3:])),
    "FTildeFunction": lambda *args: Z3num.FTilde_norm_symmetrical(*args[:-3], *map(bool, args[-3:])),
    "cF3": Z3num.cF3_norm_symmetrical,
    "F3": Z3num.F3
})
num_module_norm = [num_module_norm, "scipy"]

num_module = Z2sym.num_module[0]
num_module.update({
    "EFunction": lambda *args: Z3num.E_symmetrical(*args[:-3], *map(bool, args[-3:])),
    "FFunction": lambda *args: Z3num.F_symmetrical(*args[:-3], *map(bool, args[-3:])),
    "FTildeFunction": lambda *args: Z3num.FTilde_symmetrical(*args[:-3], *map(bool, args[-3:])),
    "cF3": Z3num.cF3_symmetrical,
    "F3": Z3num.F3
})
num_module = [num_module, "scipy"]
