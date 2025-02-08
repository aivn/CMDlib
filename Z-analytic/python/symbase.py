"""
Общие символьные конструкции.
"""
import sympy as sp


class SigmaSymbol(sp.Symbol):
    def _latex(self, printer, exp=None):
        if exp is None:
            return self.name
        else:
            exp = printer.doprint(exp)
            return f"{self.name}^{exp}"

    def _eval_power(self, expr):
        if expr.is_integer:
            return 1-sp.Mod(expr, 2) + self*sp.Mod(expr, 2)
