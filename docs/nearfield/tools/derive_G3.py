"""Canonical closed forms G_mn^(p) = P0 + P1 R + P2 asinh(U/V) + P3 asinh(V/U) (P_i polynomials in U, V) of
G_mn^(p)(U, V) = int_0^U int_0^V u^m v^n R^p dv du, R = sqrt(U^2 + V^2), for p = -1, 1, 3 and m, n = 0..3.
Verified: d2G/dUdV = U^m V^n R^p (symbolically), G -> 0 on both axes, and against 30-digit quadrature.
Writes G3.pkl (coefficient polynomials), G3_table.tex (LaTeX) and nearFieldG.m (MATLAB)."""
import pickle
import sympy as sp
import mpmath as mp

mp.mp.dps = 30
u, v, U, V = sp.symbols('u v U V', positive=True)
Rs, LU, LV = sp.symbols('R L_U L_V', positive=True)        # R, asinh(U/V), asinh(V/U)
Ru = sp.sqrt(u**2 + V**2)

def outer_integral(expr):
    total = 0
    for term in sp.Add.make_args(sp.expand(expr)):
        coeff, part = term.as_independent(u, as_Add=False)
        if part.has(sp.asinh):
            k_part = sp.simplify(part / sp.asinh(V / u))
            k = sp.degree(k_part, u) if k_part != 1 else 0
            assert sp.simplify(k_part - u**k) == 0, part
            tail = sp.integrate(u**k / Ru, (u, 0, U))
            total += coeff * (U**(k + 1) / (k + 1) * sp.asinh(V / U) + V / (k + 1) * tail)
        else:
            total += coeff * sp.integrate(part, (u, 0, U))
    return total

def canonical(expr):
    R = sp.sqrt(U**2 + V**2)
    e = expr.rewrite(sp.asinh)
    e = e.subs({sp.asinh(U / V): LU, sp.asinh(V / U): LV})
    e = e.subs(R, Rs)
    e = sp.expand(e)
    # reduce powers of R with R^2 = U^2 + V^2 and clear 1/R
    num, den = sp.fraction(sp.together(e))
    num = sp.expand(num)
    poly = sp.Poly(num, Rs)
    reduced = 0
    for (j,), c in poly.terms():
        reduced += c * (U**2 + V**2)**(j // 2) * Rs**(j % 2)
    den_poly = sp.Poly(sp.expand(den), Rs)
    assert den_poly.degree() <= 1, den
    if den_poly.degree() == 1:
        # den = d1 R + d0 with d0 = 0 expected: multiply by R / (U^2 + V^2)
        d1 = den_poly.coeff_monomial(Rs); d0 = den_poly.coeff_monomial(1)
        assert d0 == 0, den
        reduced = sp.expand(reduced * Rs)
        poly = sp.Poly(reduced, Rs)
        reduced = 0
        for (j,), c in poly.terms():
            reduced += c * (U**2 + V**2)**(j // 2) * Rs**(j % 2)
        reduced = sp.expand(reduced) / (d1 * (U**2 + V**2))
    else:
        reduced = reduced / den
    reduced = sp.expand(reduced)
    parts = {}
    rem = reduced
    for name, sym in (('P2', LU), ('P3', LV)):
        c = sp.expand(rem).coeff(sym)
        parts[name] = sp.simplify(c)
        rem = sp.expand(rem - c * sym)
    c = sp.expand(rem).coeff(Rs)
    parts['P1'] = sp.simplify(c)
    parts['P0'] = sp.simplify(sp.expand(rem - c * Rs))
    for k, val in parts.items():
        assert not val.has(LU) and not val.has(LV) and not val.has(Rs), (k, val)
        assert sp.simplify(val - sp.cancel(val)) == 0
    return parts

def rebuild(parts):
    R = sp.sqrt(U**2 + V**2)
    return parts['P0'] + parts['P1'] * R + parts['P2'] * sp.asinh(U / V) + parts['P3'] * sp.asinh(V / U)

table = {}
for p in (-1, 1, 3):
    for m in range(4):
        for n in range(4):
            inner = sp.integrate(v**n * (u**2 + v**2)**sp.Rational(p, 2), (v, 0, V))
            inner = inner.rewrite(sp.asinh) if inner.has(sp.log) else inner
            parts = canonical(outer_integral(u**m * inner))
            G = rebuild(parts)
            R = sp.sqrt(U**2 + V**2)
            check = sp.simplify(sp.diff(G, U, V) - U**m * V**n * R**p)
            polys = all(sp.cancel(parts[k]).is_polynomial(U, V) for k in parts)
            Gf = sp.lambdify((U, V), G, 'mpmath')
            err = max(abs(Gf(mp.mpf(x), mp.mpf(y)) - mp.quad(lambda s, t: s**m * t**n * (s**2 + t**2)**(mp.mpf(p) / 2), [0, x], [0, y]))
                      for (x, y) in [('0.7', '0.3'), ('0.05', '1.3'), ('2', '2')])
            lim = (sp.limit(G, U, 0), sp.limit(G, V, 0))
            table[(p, m, n)] = parts
            print(f'p={p:2d} m={m} n={n}: check={check} polynomial={polys} limits={lim} err={mp.nstr(err, 2)}', flush=True)

pickle.dump({k: {kk: sp.srepr(vv) for kk, vv in d.items()} for k, d in table.items()}, open('G3.pkl', 'wb'))

# LaTeX table (one block per p).
with open('G3_table.tex', 'w') as f:
    for p in (-1, 1, 3):
        for m in range(4):
            for n in range(4):
                d = table[(p, m, n)]
                terms = []
                for key, tail in (('P0', ''), ('P1', r'\,R'), ('P2', r'\,\operatorname{arsinh}\frac{U}{V}'), ('P3', r'\,\operatorname{arsinh}\frac{V}{U}')):
                    c = sp.factor(d[key])
                    if c != 0:
                        terms.append(r'\left(' + sp.latex(c) + r'\right)' + tail)
                f.write(r'\GLine{%d}{%d}{%d}{%s}' % (p, m, n, ' + '.join(terms) if terms else '0') + '\n')

# MATLAB: G = nearFieldG(p, m, n, U, V), elementwise, zero on the axes.
with open('nearFieldG.m', 'w') as f:
    f.write('function G = nearFieldG(p, m, n, U, V)\n')
    f.write('% nearFieldG: Closed form of G_mn^(p)(U, V) = int_0^U int_0^V u^m v^n (u^2 + v^2)^(p/2) dv du for U, V >= 0\n')
    f.write('%   (p = -1, 1, 3; m, n = 0..3), as P0 + P1 R + P2 asinh(U/V) + P3 asinh(V/U) with polynomial P_i\n')
    f.write('%   (generated by derive_G3.py, see the near-field note, Appendix A). G vanishes on both axes.\n\n')
    f.write('    R = sqrt(U.^2 + V.^2);\n')
    f.write('    key = sprintf(\'%d_%d_%d\', p, m, n);\n')
    f.write('    switch key\n')
    for (p, m, n), d in table.items():
        f.write(f"        case '{p}_{m}_{n}'\n")
        code = {k: sp.octave_code(sp.expand(d[k])) for k in d}
        f.write(f"            G = ({code['P0']}) + ({code['P1']}).*R + ({code['P2']}).*asinh(U./V) + ({code['P3']}).*asinh(V./U);\n")
    f.write('        otherwise\n            error(\'nearFieldG: no closed form for %s.\', key);\n    end\n\n')
    f.write('    G(U == 0 | V == 0) = 0;\nend\n')
print('written G3.pkl, G3_table.tex, nearFieldG.m')
