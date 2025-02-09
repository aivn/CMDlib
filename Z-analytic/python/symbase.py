"""
Общие символьные конструкции.
"""
import sympy as sp


class SigmaSymbol(sp.Symbol):
    """sigma**2 = 1"""
    def _latex(self, printer, exp=None):
        if exp is None:
            return self.name
        else:
            exp = printer.doprint(exp)
            return f"{self.name}^{exp}"

    def _eval_power(self, expr):
        if expr.is_integer:
            return 1-sp.Mod(expr, 2) + self*sp.Mod(expr, 2)


def eval_sum_direct(expr, recursive=1):
    """Рекурсивное вычисление sp.Sum."""
    if recursive < 1:
        return expr
    
    def func(expr, limits):
        if isinstance(limits[1], sp.Symbol) or isinstance(limits[2], sp.Symbol):
            return sp.Sum(expr, limits)
        return sp.concrete.summations.eval_sum_direct(expr, limits)
    
    expr = expr.replace(sp.Sum, func)
    return eval_sum_direct(expr, recursive-1)
