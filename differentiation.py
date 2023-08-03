"""package for numerical derivation"""
import numpy as np


def pascal(n=1):
    last = [1, 1]
    while len(last) - 1 < n:
        next = [1]
        for i in range(len(last) - 1):
            next.append(last[i] + last[i + 1])
        next.append(last[-1])
        last = next
    return last


def dern(fun, n, x, epsilon=5e-4):
    ps = pascal(n)
    evalpoints_n = -np.arange(len(ps)) * 2.
    evalpoints_n -= np.mean(evalpoints_n)
    results = []
    while len(results) < 1:  # or np.abs(results[-1] - results[-2]) >= 1e-7:
        evalpoints = x + epsilon * evalpoints_n
        alternating_signs = [(-1) ** i for i, _ in enumerate(ps)]
        top = sum(fun(evalpoint) * alternating_sign * p
                  for evalpoint, alternating_sign, p in zip(evalpoints, alternating_signs, ps))
        bot = (2 * epsilon) ** n
        results.append(top / bot)
        epsilon /= 2.
    return results[-1]

if __name__ == '__main__':
    print(dern(np.sin,1,np.pi/4))