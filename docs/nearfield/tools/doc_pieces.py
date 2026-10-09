"""Pieces and moments of the near-field note: the polynomial pieces of the weights in u per offset (Lemma pieces), the
touching moments as corner sums of G (Theorem touching), and the symbolic proofs of the coincident identities with the
src/thesis formulas (Proposition coincident). The part above '# (3)' is shared by make_tex.py and exact_values.py."""
import pickle
import sympy as sp

a, b = sp.symbols('a b', positive=True)
U, V = sp.symbols('U V', positive=True)
nu, Up = sp.symbols('nu u', real=True)
G3 = {k: {kk: sp.S(eval(vv, sp.__dict__)) for kk, vv in d.items()} for k, d in pickle.load(open('G3.pkl', 'rb')).items()}
def Gfun(p, m, n):
    d = G3[(p, m, n)]
    R = sp.sqrt(U**2 + V**2)
    return d['P0'] + d['P1'] * R + d['P2'] * sp.asinh(U / V) + d['P3'] * sp.asinh(V / U)

W = {
    (0, 0): {'-': 1 + nu, '+': 1 - nu},
    (1, 0): {'-': (1 + nu)**2 / 2, '+': (1 - nu**2) / 2},
    (0, 1): {'-': (1 - nu**2) / 2, '+': (1 - nu)**2 / 2},
    (1, 1): {'-': (1 + nu)**2 * (2 - nu) / 6, '+': sp.Rational(1, 3) - nu / 2 + nu**3 / 6},
}
pieces = {'-': (-1, 0), '+': (0, 1)}

def piece(wkey, side, d, h):
    n0, n1 = pieces[side]
    ends = [h * (n0 - d), h * (n1 - d)]
    sigma = 1 if min(ends) >= 0 else -1
    lo, hi = sorted([sigma * e for e in ends])
    poly = sp.expand(W[wkey][side].subs(nu, sigma * Up / h + d))
    return sigma, lo, hi, poly

def moment(p, wx, wy, dx, dy):
    total = 0
    for sx in ('-', '+'):
        for sy in ('-', '+'):
            _, u0, u1, px = piece(wx, sx, dx, a)
            _, v0, v1, py = piece(wy, sy, dy, b)
            Px = sp.Poly(px, Up); Py = sp.Poly(py.subs(Up, sp.Symbol('Vp')), sp.Symbol('Vp'))
            for (m,), cm in Px.terms():
                for (n,), cn in Py.terms():
                    g = Gfun(p, m, n)
                    val = 0
                    for (x, s1) in ((u1, 1), (u0, -1)):
                        for (y, s2) in ((v1, 1), (v0, -1)):
                            if x != 0 and y != 0:
                                val += s1 * s2 * g.subs({U: x, V: y})
                    total += cm * cn * val
    return total / (a * b)

# (3) coincident identities with the src / thesis formulas (symbolic, with a = 1 by homogeneity, b > 0 free)
L = lambda x, y: sp.log((sp.sqrt(x**2 + y**2) + x) / y)
s_ = sp.sqrt(a**2 + b**2)
type3 = 2 * (a**3 - a**2 * s_ + b**3 - b**2 * s_ + 3 * a * b**2 * L(a, b) + 3 * a**2 * b * L(b, a)) / (3 * a**2 * b**2)
type1 = (12 * a**5 - 12 * a**4 * s_ + a**2 * b**2 * (30 * b - 29 * s_) + 2 * b**4 * (b - s_) + 15 * a**3 * b * (4 * b * L(a, b) + 3 * a * L(b, a))) / (90 * a**4 * b**2)
type2 = (18 * a**5 - 2 * b**5 - 18 * a**4 * s_ - a**2 * b**2 * s_ + 2 * b**4 * s_ + 30 * a**3 * b**2 * L(a, b) + 45 * a**4 * b * L(b, a)) / (90 * a**4 * b**2)
I1one = (4 * a**5 + 4 * b**5 - s_ * (4 * a**4 - 12 * a**2 * b**2 + 4 * b**4) + 10 * a * b**4 * L(a, b) + 10 * a**4 * b * L(b, a)) / 60 / (a * b)**2
I1st = (36 * a**7 + 126 * a**2 * b**5 + 8 * b**7 - s_ * (36 * a**6 - 165 * a**4 * b**2 + 122 * a**2 * b**4 + 8 * b**6) + 210 * a**3 * b**4 * L(a, b) + 105 * a**6 * b * L(b, a)) / (3780 * a**2) / (a * b)**2
checks = [('M00^(-1)(0,0) - type3', moment(-1, (0, 0), (0, 0), 0, 0) - type3),
          ('M11^(-1)(0,0) - type1', moment(-1, (1, 1), (0, 0), 0, 0) - type1),
          ('M10^(-1)(0,0) - M11^(-1)(0,0) - type2', moment(-1, (1, 0), (0, 0), 0, 0) - moment(-1, (1, 1), (0, 0), 0, 0) - type2),
          ('M00^(1)(0,0) - thesis I1[1]', moment(1, (0, 0), (0, 0), 0, 0) - I1one),
          ('M11^(1)(0,0) - thesis I1[lambda lambda]', moment(1, (1, 1), (0, 0), 0, 0) - I1st)]
with open('identities.txt', 'w') as f:
    for name, diff in checks:
        e = sp.simplify(sp.expand(sp.expand_log(diff.rewrite(sp.log), force=True)).subs(a, 1))
        numeric = sp.N(diff.subs({a: sp.Rational(7, 10), b: sp.Rational(3, 10)}), 40)
        line = f'{name}: simplified (a = 1) = {e}; at a = 0.7, b = 0.3: {sp.N(numeric, 3)}'
        print(line, flush=True); f.write(line + '\n')

print('done')
