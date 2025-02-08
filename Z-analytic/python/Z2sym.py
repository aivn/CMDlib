"""
Символьный расчет двухчастичной статистической суммы.
"""
import sympy as sp
import Z2num

SIGMAS = [-1, 1]

h, l = sp.Symbol("h"), sp.Symbol("\\lambda", positive=True)
m, eta = sp.Symbol("\\left<m\\right>"), sp.Symbol("\\left<\eta\\right>", positive=True)


num_module_norm = {
    "cF2_symmetrical": Z2num.cF2_norm_symmetrical,
    "E_symmetrical": lambda *args: Z2num.E_norm_symmetrical(*args[:2], *map(bool, args[2:])),
    "F2": Z2num.dawson
}
num_module_norm = [num_module_norm, "scipy"]

num_module = {
    "cF2_symmetrical": Z2num.cF2_symmetrical,
    "E_symmetrical": lambda *args: Z2num.E_symmetrical(*args[:2], *map(bool, args[2:])),
    "F2": Z2num.dawson
}
num_module = [num_module, "scipy"]

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
        

class E_symmetrical(sp.Function):
    """
    args[0], args[1] -- h, l
    args[2]: {0, 1} -- есть ли sigma_1 в множителе внутри суммы
    args[3]: {0, 1} -- есть ли sigma_2 в множителе внутри суммы
    EFunction(h, l, 1, 1).diff(h) -> 2*EFunction(h, l, 1, 0)
    """

    def _latex(self, printer, exp=None):
        args = printer.doprint(self.args[0]), printer.doprint(self.args[1])
        args = f"\\left({args[0]}, {args[1]}\\right)"

        indexes = ""
        for i, item in enumerate(self.args[2:]):
            if item:
                indexes += str(i+1)
        name = "{\\cal E}_{"+indexes+"}"

        if exp is None:
            return name+args
        else:
            exp = printer.doprint(exp)
            return f"{name}^{exp}"+args

    def fdiff(self, argindex=1):
        if argindex == 1:
            args = self.args[0], self.args[1]
            term1 = E_symmetrical(
                *args,
                int(not (self.args[2] and True)),
                self.args[3]
            )
            term2 = E_symmetrical(
                *args,
                self.args[2],
                int(not (self.args[3] and True))
            )
            return term1 + term2

        if argindex == 2:
            args = self.args[0], self.args[1]
            return E_symmetrical(
                *args,
                int(not (self.args[2] and True)),
                int(not (self.args[3] and True))
            )

    @classmethod
    def eval(cls, *args):
        if (not args[2]) and args[3]:
            return E_symmetrical(args[0], args[1], 1, 0)

    def _eval_rewrite(self, rule, args, **hints):
        h, l, arg_s1, arg_s2 = args
        if rule == sp.exp:
            result = 0
            for s1 in SIGMAS:
                for s2 in SIGMAS:
                    term = sp.exp((s1+s2)*h+s1*s2*l)
                    if arg_s1:
                        term *= s1
                    if arg_s2:
                        term *= s2
                    result += term
            return result

class cF2_symmetrical(sp.Function):
    def _latex(self, printer, exp=None):
        args = printer.doprint(self.args[0]), printer.doprint(self.args[1])
        if exp is None:
            return "{\\cal F}^{(2)}"+f"\\left({args[0]},{args[1]}\\right)"
        else:
            exp = printer.doprint(exp)
            return f"\\left({self._latex(printer)}\\right)^{exp}"

    def fdiff(self, argindex=1):
        h, l = self.args
        if argindex == 1:
            return -h/l*cF2_symmetrical(h, l) + 1/sp.sqrt(2*l)*E_symmetrical(h, l, 1, 1)

        if argindex == 2:
            term1 = (h**2/(2*l**2)-1)*cF2_symmetrical(h, l)
            term2 = -h/(2*l*sp.sqrt(2*l))*E_symmetrical(h, l, 1, 1)
            term3 = 1/sp.sqrt(2*l)*E_symmetrical(h, l, 1, 0)
            return term1 + term2 + term3

    def _eval_rewrite(self, rule, args, **hints):
        h, l = args
        if rule == sp.exp:
            result = 0
            for s1 in SIGMAS:
                for s2 in SIGMAS:
                    term = sp.exp((s1+s2)*h+s1*s2*l)
                    term *= F2((h+(s1+s2)*l)/sp.sqrt(2*l))
                    term *= s1*s2
                    result += term
            return result