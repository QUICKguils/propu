import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from propu.util import isatmosphere
from propu.util.constant import ur


def main():
    print("Solving exo 2.3: MQ-1B Predator")

    # Statement data
    power = (115*ur.hp).to(ur.W)
    rho, _, _, _ = isatmosphere.get_state(0*ur.m)
    diameter = (6*ur.ft).to(ur.m)

    # Resolution
    Ap = np.pi * diameter**2 / 4

    def f(eta, v_inf):
        return eta - v_inf / (v_inf/2 + np.sqrt(v_inf**2/4 + eta*power/(2*Ap*rho*v_inf)))

    guess = 0.5
    v_range = np.arange(2, 102, 2) * ur("m/s")
    eta = np.zeros(v_range.shape)
    for (i, v) in enumerate(v_range):
        sol = optimize.root(f, guess, args=v)
        eta[i] = sol.x

    fig, ax = plt.subplots()
    ax.plot(v_range, eta)
    fig.show()

    return (v_range, eta)
