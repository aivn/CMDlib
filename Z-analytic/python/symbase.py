"""
Общие символьные конструкции.
"""
import sympy as sp
from sympy.printing.c import C99CodePrinter


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


def collect_sigmas(expr, terms, sigmas):
    """Разделяет выражение по заданным terms, после чего дополнительно разделяет по степеням заданных sigmas."""
    dummy_subs = dict(zip(terms, [sp.Dummy() for _ in range(len(terms))]))
    expr = expr.expand().subs(dummy_subs)

    result = sp.collect(expr, dummy_subs.values(), evaluate=False)
    for dummy in dummy_subs.values():
        if dummy in result:
            result[dummy] = sp.Poly(result[dummy], sigmas).as_dict()
    for key, dummy in dummy_subs.items():
        if dummy in result:
            result[key] = result.pop(dummy)

    return result


def LHopital(expr, variable, N=1):
    """Дифференцирует числитель и знаменатель у выражения по отдельности."""
    if N == 0:
        return expr
    n, d = sp.fraction(expr)
    return LHopital((n.diff(variable)/d.diff(variable)).factor(), variable, N-1)


class CPrinter(C99CodePrinter):
    """Класс для генерации C-кода."""
    def _print_Integer(self, expr):
        return str(expr) + "."

    def _print_Rational(self, expr):
        s = str(expr.p/expr.q)
        if len(s) < 14:
            return s
        else:
            return super()._print_Rational(expr)

    def _print_Pow(self, expr):
        arg, power = expr.args
        if power.is_integer:
            return "("+"*".join([self._print(arg)]*int(power.evalf()))+")"

        return super()._print_Pow(expr)


def get_ccode(expr):
    """Генерация C-кода."""
    return CPrinter().doprint(expr)


def optimize_pow(expr, names):
    """Оптимизация pow."""
    def _optimize_pow(*args):
        symbol = args[0]
        if isinstance(symbol, sp.Symbol):
            name_symbol = names.get(symbol, symbol.name)
            power = args[1]

            if power.is_integer:
                return sp.Symbol(f"{name_symbol}")**power
            if (2*power).is_integer:
                if power < 0:
                    power = int(power + 1/2)
                    return 1/sp.Symbol(f"sqrt_{name_symbol}") * sp.Symbol(f"{name_symbol}")**(power)
                if power > 0:
                    power = int(power - 1/2)
                    return sp.Symbol(f"sqrt_{name_symbol}") * symbol**(power)

    expr = expr.replace(sp.Pow, _optimize_pow)
    return expr