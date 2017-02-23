#!/usr/bin/env python3
# coding:utf-8

import sys
import itertools
from math import sqrt
from scipy.stats import norm
import numpy as np


class Combination(object):
    def __init__(self):
        self.combination = None
        self.wavelength = None
        self.k_ion = None
        self.sigma = list()
        self.p = list()

    def __str__(self):
        i = self.combination[0]
        j = self.combination[1]
        k = self.combination[2]
        message = (
            '{0} {1} {2} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7:.3f}').format(
                i, j, k, self.wavelength, self.k_ion, self.sigma[0],
                self.sigma[1], self.sigma[2])
        return message

    def set(self, i, j, k):
        self.combination = (i, j, k)
        wavelength = c / (i * f1 + j * f2 + k * f3)
        beta_ijk = f1 * f1 * (i / f1 + j / f2 + k / f3) / (
            i * f1 + j * f2 + k * f3)
        beta_lmn = l + m * f1 * f1 / (f2 * f2) + n * f1 * f1 / (f3 * f3)
        k_ion = (beta_ijk + beta_lmn) / wavelength
        self.wavelength = wavelength
        self.k_ion = k_ion
        # set sigma, p
        for sigma_p, sigma_fi in zip([0.3, 0.6, 3], [0.01, 0.01, 0.01]):
            sigma = 2 * sqrt((l * l + m * m + n * n) * sigma_p * sigma_p / (
                wavelength * wavelength) + (i * i + j * j + k * k) * sigma_fi *
                             sigma_fi)

            p = (norm.cdf(0.5, 0, sigma) - norm.cdf(-0.5, 0, sigma)) * 100
            self.sigma.append(sigma)
            self.p.append(p)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('\nArgs:')
        print('\tsystem')
        print('\tthreshold')
        sys.exit()

    c = 299792458

    if sys.argv[1] == 'b':
        f1 = 1561.098 * (10**6)
        f2 = 1207.140 * (10**6)
        f3 = 1268.520 * (10**6)
    elif sys.argv[1] == 'g':
        f1 = 1575.42 * (10**6)
        f2 = 1227.60 * (10**6)
        f3 = 1176.45 * (10**6)
    l = 1 / 3.0
    m = 1 / 3.0
    n = 1 / 3.0
    threshold = float(sys.argv[2])
    combinations = list()
    for i in range(-10, 11):
        for j in range(-10, 11):
            for k in range(-10, 11):
                if i == 0 and j == 0 and k == 0:
                    continue

                if i * f1 + j * f2 + k * f3 == 0:
                    continue

                if abs(i + j + k) > 2:
                    continue

                combination = Combination()
                combination.set(i, j, k)
                if combination.wavelength > 4 and combination.sigma[
                        2] < threshold:
                    combinations.append(combination)
                    print(combination)

    print('\nGroup:')
    # start search group
    for combination in itertools.combinations(combinations, 3):
        min_p = list()
        for i in range(3):
            temp = list()
            for item in combination:
                temp.append(item.p[i])
            min_p.append(min(temp))

        if min_p[0] > 99.9 and min_p[1] > 99.9 or min_p[2] > 94:
            group = [item.combination for item in combination]
            matrix = np.array(group)
            if abs(abs(np.linalg.det(matrix)) - 1) < 0.001:
                print(* [x for item in group for x in item],
                      * ['{:.2f}%'.format(p) for p in min_p])
