from itertools import product
from collections import defaultdict
from tqdm import tqdm
from sage.symbolic import operators
from sage.symbolic.expression_conversions import FakeExpression

class Monomial:
    """Hack to make sure a monomial is treated as an addition in the AST, not a
    multiplciation."""
    def __init__(self, expr):
        self.expr = expr

    def variables(self):
        return self.expr.variables()

    def is_polynomial(self, var):
        return self.expr.is_polynomial(var)

    def expand(self):
        return Monomial(self.expr.expand())

    def operands(self):
        return [self.expr]

    def __neg__(self):
        return Monomial(-self.expr)


# https://ask.sagemath.org/question/7890/solve-for-all-values-of-a-variable/
def generate_polynomial_constraints(lhs, rhs, coeff_vars):
    # get the non-coefficient variables
    vv = sorted(set(lhs.variables()).union(set(rhs.variables())).difference(coeff_vars))

    # make sure they're polynomials
    assert all(lhs.is_polynomial(v) and rhs.is_polynomial(v) for v in vv)
    d = defaultdict(Integer)
    # loop over each term in the equation
    for term in (lhs.expand().operands())+(-(rhs.expand())).operands():
        # get the degrees of the variables
        vsig = tuple(term.degree(v) for v in vv)
        # find the non-variable part
        d[vsig] += term/(prod(v**vs for v,vs in zip(vv, vsig)))
    # get solutions with free variables
    return d.values()

def generate_rational_constraints(eq, coeff_vars):
    l_num = eq.lhs().numerator().expand()
    r_num = eq.rhs().numerator().expand()
    l_den = eq.lhs().denominator().expand()
    r_den = eq.rhs().denominator().expand()

    if l_num.operator() == operators.mul_vararg:
        l_num = Monomial(l_num)

    if r_num.operator() == operators.mul_vararg:
        r_num = Monomial(r_num)

    if l_den.operator() == operators.mul_vararg:
        l_den = Monomial(l_den)

    if r_den.operator() == operators.mul_vararg:
        r_den = Monomial(r_den)

    print(l_num, r_num)

    return [
        *generate_polynomial_constraints(l_num, r_num, coeff_vars),
        *generate_polynomial_constraints(l_den, r_den, coeff_vars)]
# x, y, z = var("x y z")
n = var("n")
e1(x, y, z) = x + y + z
e2(x, y, z) = x*y + x*z + y*z
e3(x, y, z) = x*y*z

# t, u, v = var("t u v")

def array_to_polynomial(arr):
    return sum(
        sum(
            sum(
                arr[i][j][k] * t ** i * u ** j * v ** k
                for k in range(len(arr[i][j]))
                )
            for j in range(len(arr[i]))
            )
        for i in range(len(arr))
    )

def check_ratl(f):
    return bool(
        f.subs(z == 1) == x*y
        and
        f.derivative(z).subs(z == 1) == x*y/(x+y-x*y)
    )

def solve_soln(f, coeff_vars):
    A = generate_rational_constraints((f.subs(z == 1) == x * y).simplify_full(), coeff_vars)
    B = generate_rational_constraints(f.derivative(z).subs(z == 1).simplify_full() == (x * y / (x + y - x * y)).simplify_full(), coeff_vars)

    print(A)
    print(B)

    return solve(B, *coeff_vars)


def iterate(degree, lower_coeff, higher_coeff):
    for data in product(range(lower_coeff, higher_coeff + 1), repeat=3 ** degree):
        pol = sum(
            data[i * degree ** 2 + j * degree + k] * t ** i * u ** j * v ** k
            for (i, j, k) in product(range(degree), repeat=3)
        )

        yield pol

def symmetrize(f):
    return f.subs(t == e1, u == e2, v == e3)

# DEGREE = 2
# BOUNDS = (-2, 2)
#
# for p in iterate(DEGREE, *BOUNDS):
#     for q in iterate(DEGREE, *BOUNDS):
#         if q != 0:
#             pol = p / q * x*y*z
#             try:
#                 if check_ratl(symmetrize(pol)):
#                     print("suc", symmetrize(pol).simplify_full())
#                     print()
#             except:
#                 print("err", symmetrize(pol))
#



# g = x*y*z * (-1+x+y+z - x*y*z) / (x*y+x*z+y*z-2*x*y*z)
# g2 = x*y*z * (2 * (x+y+z) - (x*y+x*z+y*z) - 2)/(x+y+z-1-x*y*z)
#
# exit()
#
# print(check_ratl(g2))

# a, b, c, d, h, i, j, k = var("a b c d h i j k")

def with_coefficients(a, b, c, d, h, i, j, k):
    # l = 1
    # m = 1
    # n = 1
    # r = 1
    # s = 1
    # t = 1
    # return e3 * (a * e1 ** l + b * e2 ** m + c * e3 ** n + d) / (h * e1 ** r + i * e2 ** s + j * e3 ** t + k)
    return e3 * (a * e1 + b * e2 + c * e3 + d) / (h * e1 + i * e2 + j * e3 + k)

# for tup in [
#     (-2, 1, 0, 2, -1, 0, 1, 1),
#     (-1, 0, 1, 1, 0, -1, 2, 0),
#     (1, 0, -1, -1, 0, 1, -2, 0),
#     (2, -1, 0, -2, 1, 0, -1, -1),
#     ]:
#     print(with_coefficients(*tup))
#
# exit()

# BOUNDS = (-2, 4)
# with tqdm(total=int((BOUNDS[1] - BOUNDS[0]) ** 8)) as pbar:
#     for tup in product(range(*BOUNDS), repeat=8):
#         try:
#             f = with_coefficients(*tup)
#             if check_ratl(f):
#                 print(tup)
#                 print(f)
#                 print()
#         except:
#             pass
#         finally:
#             pbar.update(int(1))
#
# exit()

# guess = guess.subs(b == 1 - a, c == a - 2, d == -a, h == a - 1, i == 2 - a, j == a - 3, k == -a + 1)

# print(solve_soln(guess, [a, b, c, d, h, i, j, k]))
# exit()

g = e3 * (n * e1 - (n+1) * e2 + (n+2) * e3 - n) / ((n+1) * e1 - (n+2) * e2 + (n+3) * e3 - (n+1))
# print(check_ratl(g))

def combine(i, j, k, l, m):
    return l * g.subs(n == i) + m * g.subs(n == j) + (1-l-m) * g.subs(n == k)


i, j, k, a, b = var("i j k a b")
g2 = combine(i, j, k, a, b)



# print(check_ratl(g.simplify_full()))

def exponentialize(f):
    return f.subs(x == e ** x, y == e ** y, z == e ** z)

def check_full(f):
    return bool(
        f.subs(z == 0) == e^(x+y)
        and
        f.derivative(z).subs(z == 0) == e^(x+y)/(e^x+e^y-e^(x+y))
    )

f = exponentialize(g)
f2 = exponentialize(g2)
# print(check_full(f))
# a = -1
# f2 = f.subs(x == a * x, y == a * y, z == a * z)
# print(f2.derivative(z).subs(z == 0).simplify_full())
# exit()
#
# f = e ** (x+y+z) * (e ** x+e ** y+e ** z-1-e ** (x+y+z)) / (e ** (x+y)+e ** (x+z)+e ** (y+z)-2 * e ** (x+y+z))

# print(check_ratl(g2))

# print(f2.subs(z == 0).simplify_full())

def diff(f, m, n):
    f = f.derivative(x)
    for _ in range(m):
        f = f.derivative(y)
    for _ in range(n):
        f = f.derivative(z)

    return f.subs(y == 0, z == 0).simplify_full()

def extract(f, l, m, n):
    t = f.taylor(x, 0, l).taylor(y, 0, m).taylor(z, 0, n)
    return ((t.coefficient(x, l) * factorial(l)).coefficient(y, m) * factorial(m)).coefficient(z, n) * factorial(n)

# print(extract(f, 2, 2, 2))
# print(solve([
#     extract(f2, 2, 2, 2) == 122,
#     extract(f2, 2, 2, 3) == 898
# ], a, b, c))
print(diff(f2, 2, 2))
print(diff(f2, 2, 3))
print(diff(f2, 2, 4))
print(diff(f2, 3, 3))
# print(diff(f2, 2, 2))
# print(diff(f2, 4, 3))
#
