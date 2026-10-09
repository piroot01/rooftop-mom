"""Exact values of the seven touching moments of R^p, p = -1, 1, 3, for the tests of
rooftopMoM.assembler.singularIntegrals.generatePowerMoments, written to exactMoments.csv next to this file.

Each moment is the two-dimensional integral of the near-field note (docs/nearfield, Two-dimensional representation),

    M_ab(d) = int_{-1}^{1} int_{-1}^{1} W_ab(nu) W_00(eta) rho_d(nu, eta)^p d eta d nu,
    rho_d = sqrt(a^2 (nu - dx)^2 + b^2 (eta - dy)^2),

and N_ab the same with the weights exchanged. It is integrated by tanh-sinh quadrature in 40-digit arithmetic on the
quadrants of [-1, 1]^2, where the weights are polynomials and the singular point (dx, dy) is a corner. The closed forms
are not used, so the values check them independently.

    uv run --no-project --with mpmath python generateExactMoments.py [output.csv]
"""
import concurrent.futures
import os
import sys

import mpmath as mp

mp.mp.dps = 40

CELLS = [('0.1', '0.1'), ('0.1', '0.025'), ('0.03', '0.11'), ('1', '3'), ('1', '0.1')]
OFFSETS = [(0, 0), (1, 0), (0, 1), (1, 1)]
# Weights along x and y of M00, M10, M01, M11, N10, N01 and N11.
KINDS = [('00', '00'), ('10', '00'), ('01', '00'), ('11', '00'), ('00', '10'), ('00', '01'), ('00', '11')]


def weight(kind, nu):
    """Correlation weight W_kind of the rooftop profiles (Weights of the rooftop moments)."""
    if kind == '00':
        return 1 - abs(nu)
    if kind == '10':
        return (1 + nu)**2 / 2 if nu < 0 else (1 - nu**2) / 2
    if kind == '01':
        return (1 - nu**2) / 2 if nu < 0 else (1 - nu)**2 / 2
    return (1 - abs(nu))**2 * (2 + abs(nu)) / 6


def moment(a, b, power, offset, kinds):
    kindX, kindY = kinds
    dx, dy = offset

    def integrand(nu, eta):
        return weight(kindX, nu) * weight(kindY, eta) * mp.sqrt(a**2 * (nu - dx)**2 + b**2 * (eta - dy)**2)**power

    return sum(mp.quad(integrand, x, y) for x in ([-1, 0], [0, 1]) for y in ([-1, 0], [0, 1]))


def row(case):
    """One line of the table: the cell, the power, the offset and the seven moments to 30 digits."""
    a, b, power, offset = case
    values = [moment(mp.mpf(a), mp.mpf(b), power, offset, kinds) for kinds in KINDS]
    return ','.join([a, b, str(power), str(offset[0]), str(offset[1])] + [mp.nstr(value, 30) for value in values])


def main(path):
    cases = [(a, b, power, offset) for a, b in CELLS for power in (-1, 1, 3) for offset in OFFSETS]

    # A moment takes seconds to tens of seconds, so the rows run in parallel.
    with concurrent.futures.ProcessPoolExecutor() as pool:
        lines = list(pool.map(row, cases))

    with open(path, 'w') as output:
        output.write('a,b,power,dx,dy,M00,M10,M01,M11,N10,N01,N11\n')
        output.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                             'exactMoments.csv'))
