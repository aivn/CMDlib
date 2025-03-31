"""
Символьный расчет двухчастичной статистической суммы.
"""
import sympy as sp
import Z2num
import symbase

SIGMAS = [-1, 1]

h, l = sp.symbols("h \\lambda", positive=True)
hi, hj = sp.symbols("h_i h_j", positive=True)
s1, s2 = symbase.SigmaSymbol("\\sigma_1"), symbase.SigmaSymbol("\\sigma_2")

m, eta = sp.symbols("\\left<m\\right> \\left<\eta\\right>")
m_eta, upsilon = sp.symbols("{\\left<m\\eta\\right>} \\Upsilon")
eta2 = sp.Symbol("{\\left<\\eta^2\\right>}")
mhj = sp.Symbol("\\left<m_{h}^2\\right>")


# ===================================== Модули для sp.lambdify =====================================
num_module_norm = {
    "cF2_symmetrical": Z2num.cF2_norm_symmetrical,
    "E_symmetrical": lambda *args: Z2num.E_norm_symmetrical(*args[:-2], *map(bool, args[-2:])),
    "cF2_not_symmetrical": Z2num.cF2_norm_not_symmetrical,
    "E_not_symmetrical": lambda *args: Z2num.E_norm_not_symmetrical(*args[:-2], *map(bool, args[-2:])),
    "F2": Z2num.dawson
}
num_module_norm = [num_module_norm, "scipy"]

num_module = {
    "cF2_symmetrical": Z2num.cF2_symmetrical,
    "E_symmetrical": lambda *args: Z2num.E_symmetrical(*args[:-2], *map(bool, args[-2:])),
    "cF2_not_symmetrical": Z2num.cF2_not_symmetrical,
    "E_not_symmetrical": lambda *args: Z2num.E_not_symmetrical(*args[:-2], *map(bool, args[-2:])),
    "F2": Z2num.dawson
}
num_module = [num_module, "scipy"]


# ===================================== F2 =====================================
class F2(sp.Function):
    def _latex(self, printer, exp=None):
        arg = printer.doprint(self.args[0])
        if exp is None:
            return "F^{(2)}"+f"\\left({arg}\\right)"
        else:
            exp = printer.doprint(exp)
            return f"\\left({self._latex(printer)}\\right)^{exp}"

    def fdiff(self, argindex=1):
        return 1-2*self.args[0]*F2(self.args[0])

    def _eval_rewrite(self, rule, args, **hints):
        x, = args
        if rule == sp.erfi:
            return sp.sqrt(sp.pi)/2*sp.exp(-x**2)*sp.erfi(x)

    @classmethod
    def eval(cls, z):
        if z.could_extract_minus_sign():
            return -cls(-z)
        if z == 0:
            return 0
        

def dawsn_to_series_oo(expr, N):
    """Асимптотика функции Доусона на бесконечности."""
    def func(*args):
        z = args[0]
        result = 0
        for n in range(N):
            result += sp.factorial2(2*n-1)/(2**(n+1)*z**(2*n+1))
        return result
    return expr.replace(F2, func)


# ===================================== Для symbase.collect_sigmas =====================================
_sigmas_expr_E_symmetrical = sp.exp(h*(s1+s2)+l*s1*s2).expand()
_sigmas_expr_cF2_symmetrical = F2((h+(s1+s2)*l)/sp.sqrt(2*l)).expand()
_sigmas_terms_symmetrical = [
    _sigmas_expr_E_symmetrical*_sigmas_expr_cF2_symmetrical,
    _sigmas_expr_E_symmetrical
]

_sigmas_expr_E_not_symmetrical = sp.exp(hi*s1+hj*s2+l*s1*s2).expand()
_sigmas_expr_cF2_not_symmetrical = F2((hj+s1*l+s2*l*hj/hi)/sp.sqrt(2*l*hj/hi)).expand()
_sigmas_terms_not_symmetrical = [
    _sigmas_expr_E_not_symmetrical*_sigmas_expr_cF2_not_symmetrical,
    _sigmas_expr_E_not_symmetrical
]

# ===================================== SigmasFunction =====================================
class Z2SigmasFunctionBase(sp.Function):
    """Базовый класс для функций с суммой по двум (!) sigma (последние два аргумента)."""

    def __init__(self, name, *args, **kwargs):
        self.name = name

    def _latex(self, printer, exp=None):
        args = list(map(printer.doprint, self.args[:-2]))
        args = f"\\left({', '.join(args)}\\right)"

        indexes = ""
        for i, item in enumerate(self.args[-2:]):
            if item:
                indexes += str(i+1)
        name = f"{self.name}"+"_{" + indexes + "}"

        if exp is None:
            return name+args
        else:
            exp = printer.doprint(exp)
            return f"{name}^{exp}"+args

    def get_sigmas_factor(self):
        factor = 1
        factor *= s1 if self.args[-2] else 1
        factor *= s2 if self.args[-1] else 1
        return factor


class E_symmetrical(Z2SigmasFunctionBase):
    def __init__(self, *args, **kwargs):
        super().__init__("{\\cal E}", *args, **kwargs)

    def fdiff(self, argindex=1):
        factor = self.get_sigmas_factor()
        expr = (factor*_sigmas_expr_E_symmetrical).diff([h, l][argindex-1])
        terms = symbase.collect_sigmas(expr, _sigmas_terms_symmetrical, [s1, s2])

        args_subs = dict(zip([h, l], self.args[:-2]))

        result = 0
        for key, item in terms[_sigmas_expr_E_symmetrical].items():
            result += item.subs(args_subs).factor()*E_symmetrical(*self.args[:-2], *key)

        return result
        
    def _eval_rewrite(self, rule, args, **hints):
        args_subs = dict(zip([h, l], self.args[:-2]))
        if rule == sp.exp:
            _sigmas_expr = self.get_sigmas_factor()*_sigmas_expr_E_symmetrical
            result = 0
            for is1 in SIGMAS:
                for is2 in SIGMAS:
                    is_subs = dict(zip([s1, s2], [is1, is2]))
                    result += _sigmas_expr.subs(args_subs).subs(is_subs)
            return result

    @classmethod
    def eval(cls, *args):
        if (not args[-2]) and args[-1]:
            return E_symmetrical(*args[:-2], 1, 0)


class E_not_symmetrical(Z2SigmasFunctionBase):
    def __init__(self, *args, **kwargs):
        super().__init__("{\\cal E}", *args, **kwargs)

    def fdiff(self, argindex=1):
        factor = self.get_sigmas_factor()
        expr = (factor*_sigmas_expr_E_not_symmetrical).diff([hi, hj, l][argindex-1])
        terms = symbase.collect_sigmas( expr, _sigmas_terms_not_symmetrical, [s1, s2])

        args_subs = dict(zip([hi, hj, l], self.args[:-2]))

        result = 0
        for key, item in terms[_sigmas_expr_E_not_symmetrical].items():
            result += item.subs(args_subs).factor()*E_not_symmetrical(*self.args[:-2], *key)

        return result
        
    def _eval_rewrite(self, rule, args, **hints):
        args_subs = dict(zip([hi, hj, l], self.args[:-2]))
        if rule == sp.exp:
            _sigmas_expr = self.get_sigmas_factor()*_sigmas_expr_E_not_symmetrical
            result = 0
            for is1 in SIGMAS:
                for is2 in SIGMAS:
                    is_subs = dict(zip([s1, s2], [is1, is2]))
                    result += _sigmas_expr.subs(args_subs).subs(is_subs)
            return result

    @classmethod
    def eval(cls, *args):
        if args[0] == args[1]:
            return E_symmetrical(args[0], *args[2:])


# ===================================== cF2 =====================================
class cF2_symmetrical(sp.Function):
    def _latex(self, printer, exp=None):
        args = list(map(printer.doprint, self.args))
        if exp is None:
            return "{\\cal F}^{(2)}"+f"\\left({', '.join(args)}\\right)"
        else:
            exp = printer.doprint(exp)
            return f"\\left({self._latex(printer)}\\right)^{exp}"

    def fdiff(self, argindex=1):
        expr = (s1*s2*_sigmas_terms_symmetrical[0]).diff([h, l][argindex-1])
        sigmas_terms = symbase.collect_sigmas( expr, _sigmas_terms_symmetrical, [s1, s2])

        args_subs = dict(zip([h, l], self.args))

        result = 0
        for key, item in sigmas_terms[_sigmas_expr_E_symmetrical].items():
            result += item.subs(args_subs).factor()*E_symmetrical(*self.args, *key)
        
        for key, item in sigmas_terms[_sigmas_expr_E_symmetrical*_sigmas_expr_cF2_symmetrical].items():
            if key != (1, 1):
                raise RuntimeError
            result += item.subs(args_subs).factor()*cF2_symmetrical(*self.args)

        return result

    def _eval_rewrite(self, rule, args, **hints):
        args_subs = dict(zip([h, l], self.args))
        if rule == sp.exp:
            _sigmas_expr = s1*s2*_sigmas_expr_E_symmetrical*_sigmas_expr_cF2_symmetrical
            result = 0
            for is1 in SIGMAS:
                for is2 in SIGMAS:
                    is_subs = dict(zip([s1, s2], [is1, is2]))
                    result += _sigmas_expr.subs(args_subs).subs(is_subs)
            return result
        

class cF2_not_symmetrical(sp.Function):
    def _latex(self, printer, exp=None):
        args = list(map(printer.doprint, self.args))
        if exp is None:
            return "{\\cal F}^{(2)}"+f"\\left({', '.join(args)}\\right)"
        else:
            exp = printer.doprint(exp)
            return f"\\left({self._latex(printer)}\\right)^{exp}"

    def fdiff(self, argindex=1):
        expr = (s1*s2*_sigmas_terms_not_symmetrical[0]).diff([hi, hj, l][argindex-1])
        sigmas_terms = symbase.collect_sigmas( expr, _sigmas_terms_not_symmetrical, [s1, s2])

        args_subs = dict(zip([hi, hj, l], self.args))

        result = 0
        for key, item in sigmas_terms[_sigmas_expr_E_not_symmetrical].items():
            result += item.subs(args_subs).factor()*E_not_symmetrical(*self.args, *key)
        
        for key, item in sigmas_terms[_sigmas_expr_E_not_symmetrical*_sigmas_expr_cF2_not_symmetrical].items():
            if key != (1, 1):
                raise RuntimeError
            result += item.subs(args_subs).factor()*cF2_not_symmetrical(*self.args)

        return result

    def _eval_rewrite(self, rule, args, **hints):
        args_subs = dict(zip([hi, hj, l], self.args))
        if rule == sp.exp:
            _sigmas_expr = s1*s2*_sigmas_expr_E_not_symmetrical*_sigmas_expr_cF2_not_symmetrical
            result = 0
            for is1 in SIGMAS:
                for is2 in SIGMAS:
                    is_subs = dict(zip([s1, s2], [is1, is2]))
                    result += _sigmas_expr.subs(args_subs).subs(is_subs)
            return result
    
    @classmethod
    def eval(cls, *args):
        if args[0] == args[1]:
            return cF2_symmetrical(args[0], args[2])
    
