"""LaTeX of the near-field note in the form of the thesis appendix: one integer prefactor 1/D, integer polynomials in
a and b, sqrt and log terms with polynomial coefficients, continuation lines aligned with \\phantom. Writes
G_table.tex (Appendix A: G_mn^(p)(a, b), m <= 3, n <= 1), explicit_table.tex (Appendix B: the touching moments) and
pieces_tabular.tex (Table of Lemma pieces). Every formula is rebuilt from its printed parts and checked to 30 digits."""
import re
import sympy as sp
exec(open('doc_pieces.py').read().split('# (3) coincident identities')[0])   # a, b, U, V, Up, Gfun, W, piece, moment

WIDTH = 38                     # characters of a display line, measured by visible()
PROBE = {a: sp.Rational(7, 10), b: sp.Rational(3, 10)}


def group(s, i):
    """The braced group starting at s[i] == '{' and the index after it."""
    depth = 0
    for j in range(i, len(s)):
        depth += {'{': 1, '}': -1}.get(s[j], 0)
        if depth == 0:
            return s[i + 1:j], j + 1


def visible(s):
    """Rough printed width of a LaTeX string in characters: a fraction is as wide as its wider part, a superscript is
    narrow, and a log or a root adds its symbol."""
    width, i = 0.0, 0
    while i < len(s):
        if s.startswith(r'\frac', i):
            top, i = group(s, i + 5)
            bottom, i = group(s, i)
            width += max(visible(top), visible(bottom)) * 0.8 + 0.5
        elif s.startswith(r'\sqrt', i):
            body, i = group(s, i + 5)
            width += visible(body) + 1.5
        elif s.startswith(r'\log', i):
            width, i = width + 3, i + 4
        elif s[i] == '^':
            body, i = group(s, i + 1) if s[i + 1] == '{' else (s[i + 1], i + 2)
            width += 0.6 * visible(body)
        elif s[i] == '\\':
            i += len(re.match(r'\\[a-zA-Z]+|\\.', s[i:]).group())
        elif s[i] in '{} _':
            i += 1
        else:
            width, i = width + (1.6 if s[i] in '+-=' else 1), i + 1
    return width


def mono(i, j):
    return ' '.join(f'{x}^{{{k}}}' if k > 1 else x for x, k in (('a', i), ('b', j)) if k > 0)


def poly_terms(P):
    """Signed terms of an integer polynomial in a, b, highest power of a first."""
    out = []
    for (i, j), c in sorted(sp.Poly(P, a, b).terms(), key=lambda t: (-t[0][0], t[0][1])):
        m = mono(i, j)
        num = '' if abs(c) == 1 and m else str(abs(c))
        out.append(('-' if c < 0 else '+') + ' ' + (num + ' ' + m).strip())
    return out


def content(Q):
    """Q = c * monomial * Q' with c the signed integer content and Q' primitive (positive leading coefficient)."""
    poly = sp.Poly(Q, a, b)
    c, prim = poly.primitive()
    i0 = min(e[0] for e in prim.monoms()); j0 = min(e[1] for e in prim.monoms())
    rest = sp.expand(prim.as_expr() / (a**i0 * b**j0))
    if sp.Poly(rest, a, b).LC() < 0:
        c, rest = -c, -rest
    return c, mono(i0, j0), rest


def factor_tex(Q, tail):
    """Signed chunk Q * tail, e.g. '- 2 a \\sqrt{a^2+b^2} (6 a^2 + b^2)'."""
    c, m, rest = content(Q)
    num = '' if abs(c) == 1 and (m or rest != 1 or tail) else str(abs(c))
    body = ' '.join(t for t in (num, m) if t)
    if rest != 1:
        terms = poly_terms(rest)
        inner = ' '.join(terms)[2:] if terms[0][0] == '+' else ' '.join(terms)
        body += (' ' if body else '') + (tail + r'\left(' + inner + r'\right)' if tail.startswith(r'\sqrt') else r'\left(' + inner + r'\right) ' + tail)
    else:
        body += (' ' if body else '') + tail
    return ('-' if c < 0 else '+') + ' ' + body.strip()


def log_tex(arg):
    """arsinh(x/y) = log((sqrt(x^2+y^2) + x)/y), the form of the thesis."""
    x, y = sp.fraction(arg)
    return r'\log\left(\frac{%s + %s}{%s}\right)' % (sp.latex(sp.sqrt(sp.expand(x**2 + y**2))), sp.latex(x), sp.latex(y))


def split(e):
    """e = (P0 + sum_i Q_i sqrt(B_i) + sum_j Q_j arsinh(A_j)) / D with integer polynomials and an integer monomial D."""
    e = sp.expand(e.replace(lambda x: x.is_Pow and x.exp.is_Rational and x.exp.q == 2,
                            lambda x: sp.sqrt(sp.factor_terms(x.base))**(2 * x.exp)))
    logs = sorted(e.atoms(sp.asinh), key=sp.default_sort_key)
    roots = sorted({x for x in e.atoms(sp.Pow) if x.exp == sp.S.Half}, key=sp.default_sort_key)
    L = sp.symbols(f'L0:{len(logs)}'); S = sp.symbols(f'S0:{len(roots)}')
    f = sp.expand(e.subs(dict(zip(logs, L))).subs(dict(zip(roots, S))))
    f = f.replace(lambda x: x.is_Pow and x.base in S, lambda x: (roots[S.index(x.base)].base)**(x.exp // 2) * x.base**(x.exp % 2))
    num, den = sp.fraction(sp.cancel(sp.together(sp.expand(f))))
    lcm, P = sp.Poly(num, a, b, *L, *S).clear_denoms()
    dc, dm = (den * lcm).as_coeff_Mul()
    g = sp.gcd(P.content(), dc)
    num, den = sp.expand(P.as_expr() / g), dc / g * dm
    parts = []
    for sym, atom, kind in [*zip(S, roots, ['sqrt'] * len(S)), *zip(L, logs, ['log'] * len(L))]:
        q = num.coeff(sym)
        assert not q.free_symbols & set(L + S), (atom, q)
        if q != 0:
            parts.append((kind, atom, q))
        num = sp.expand(num - q * sym)
    assert not num.free_symbols - {a, b}
    rebuilt = (num + sum(q * atom for _, atom, q in parts)) / den
    assert abs(sp.N((rebuilt - e).subs(PROBE), 40)) < sp.Float('1e-30'), e
    return den, num, parts


def split_chunk(c):
    """A chunk wider than a line, coefficient (polynomial) tail, as pieces that break between the polynomial terms."""
    if visible(c) <= WIDTH or r'\left(' not in c:
        return [c]
    head, rest = c.split(r'\left(', 1)
    inner, tail = rest.split(r'\right)', 1)
    terms = re.split(r' (?=[+-] )', inner)
    terms[0], terms[-1] = head + r'\bigl(' + terms[0], terms[-1] + r'\bigr)' + tail
    return terms


def display(lhs, e):
    """align* rows of 'lhs = e' in the form of the thesis appendix."""
    den, P0, parts = split(e)
    chunks = poly_terms(P0) if P0 != 0 else []
    chunks += [factor_tex(q, sp.latex(atom)) for kind, atom, q in parts if kind == 'sqrt']
    chunks += [factor_tex(q, log_tex(atom.args[0])) for kind, atom, q in parts if kind == 'log']
    open_, close = (r'\frac{1}{%s}\Bigg(' % sp.latex(den), r'\Bigg)') if den != 1 else ('', '')
    lines, line = [], ''
    for c in [piece for chunk in chunks for piece in split_chunk(chunk)]:
        if line and visible(line) + visible(c) > WIDTH:
            lines.append(line); line = ''
        line += (' ' if line else '') + c
    lines.append(line)
    lines[0] = lines[0][2:] if lines[0].startswith('+ ') else lines[0]
    rows = [lhs + r' &= ' + open_ + lines[0]]
    rows += [r'&\phantom{{}=' + open_ + '}' + l for l in lines[1:]]
    rows[-1] += close
    return (r'\\' + '\n').join(rows)


def block(entries):
    return r'\begin{align*}' + '\n' + (r'\\[1.2ex]' + '\n').join(entries) + '\n' + r'\end{align*}' + '\n'


def piece_tex(poly):
    """A weight on a piece as one fraction of factors with a positive leading term in a, e.g. (a - u)^2 (2a + u)/(6a^3)."""
    num, den = sp.fraction(sp.cancel(poly))
    c, factors = sp.factor_list(num)
    out = []
    for f, k in factors:
        t = sp.Poly(f, a, Up).terms()
        if max(t, key=lambda x: (x[0][0], -x[0][1]))[1] < 0:
            f, c = -f, c * (-1)**k
        s = ' '.join(poly_terms_in(f))
        s = s[2:] if s.startswith('+ ') else s
        bracket = len(t) > 1 and (len(factors) > 1 or k > 1)
        out.append((r'\left(%s\right)' % s if bracket else s) + (f'^{{{k}}}' if k > 1 else ''))
    body = ' '.join(out)
    c, den = sp.Rational(c) / sp.Poly(den, a).LC(), den / sp.Poly(den, a).LC()
    num_c, den_c = sp.fraction(c)
    sign = '- ' if num_c < 0 else ''
    top = (str(abs(num_c)) + ' ' if abs(num_c) != 1 else '') + body
    return sign + (r'\frac{%s}{%s}' % (top, sp.latex(den_c * den)) if den_c * den != 1 else top)


def poly_terms_in(f):
    """Signed terms of a polynomial in a and u, highest power of a first."""
    out = []
    for (i, j), c in sorted(sp.Poly(f, a, Up).terms(), key=lambda t: (-t[0][0], t[0][1])):
        m = ' '.join(f'{x}^{{{k}}}' if k > 1 else x for x, k in (('a', i), ('u', j)) if k > 0)
        num = '' if abs(c) == 1 and m else str(abs(c))
        out.append(('-' if c < 0 else '+') + ' ' + (num + ' ' + m).strip())
    return out


ab = {U: a, V: b}
with open('G_table.tex', 'w') as f:
    for p in (-1, 1, 3):
        f.write(r'\paragraph{Kernel $R^{%d}$.}' % p + '\n')
        f.write(block([display(r'G^{(%d)}_{%d%d}(a,b)' % (p, m, n), Gfun(p, m, n).subs(ab, simultaneous=True))
                       for m in range(4) for n in range(2)]))

names = {(0, 0): '00', (1, 0): '10', (0, 1): '01', (1, 1): '11'}
labels = {(0, 0): r'Coincident cell, $\offs=(0,0)$', (1, 0): r'Shared edge parallel to $y$, $\offs=(1,0)$',
          (0, 1): r'Shared edge parallel to $x$, $\offs=(0,1)$', (1, 1): r'Shared vertex, $\offs=(1,1)$'}
with open('explicit_table.tex', 'w') as f:
    for p in (-1, 1, 3):
        f.write(r'\subsection*{Kernel $R^{%d}$}' % p + '\n')
        for dd in labels:
            f.write(r'\paragraph{%s.}' % labels[dd] + '\n')
            f.write(block([display(r'M^{(%d)}_{%s}(%d,%d)' % (p, label, *dd), moment(p, wkey, (0, 0), *dd))
                           for wkey, label in names.items()]))

with open('pieces_tabular.tex', 'w') as f:
    cols = [(d, side) for d in (0, 1) for side in ('-', '+')]
    heads = [piece((0, 0), side, d, a) for d, side in cols]
    f.write(r'\begin{tabular}{c' + 'c' * len(cols) + '}\n' + r'    \toprule' + '\n')
    f.write('    $d_x$ & ' + ' & '.join(r'$%d$' % d for d, _ in cols) + r' \\' + '\n')
    f.write('    piece & ' + ' & '.join(r'$I_{%s}$, $\sigma=%s$' % (side, '+1' if h[0] > 0 else '-1')
                                   for (d, side), h in zip(cols, heads)) + r' \\' + '\n')
    f.write('    $u$ & ' + ' & '.join(r'$[%s,\,%s]$' % (sp.latex(h[1]), sp.latex(h[2])) for h in heads)
            + r' \\ \midrule' + '\n')
    for wkey in W:
        f.write(r'    $W_{%d%d}$ & ' % wkey + ' & '.join(r'$\displaystyle %s$' % piece_tex(piece(wkey, side, d, a)[3])
                                                      for d, side in cols) + r' \\[2ex]' + '\n')
    f.write(r'    \bottomrule' + '\n' + r'\end{tabular}' + '\n')
print('written G_table.tex, explicit_table.tex, pieces_tabular.tex')
