import propu.constant as c
from propu.iteralg import IterTable


def main():
    print("Solving exo 3.2: Ten-stage compressor")

    # Statement data
    mdot = 164 * c.uconv("lb/s", "kg/s")
    pi_stage = 1.3
    n_stage = 10
    eta_p = 0.9
    # p0_in = 57_452  # [Pa]
    T0_in = 266.42  # [K]

    # Instantiate a table printer
    variables = ("Iter", "Power", "eta_s_tot", "eta_s_stage", "gamma")
    units = ("", "(MW)", "", "(%)", "(%)")
    table = IterTable(variables, units)

    # Total pressure ratio over the ten stages
    pi_tot = pi_stage**n_stage

    # Iterations on gamma value.
    # cp needs T0_out, which needs gamma, which in turns need cp.
    n_iter = 5
    g = c.gamma_air
    for iter in range(1, n_iter + 1):
        T0_out = T0_in * pi_tot ** ((g - 1) / (eta_p * g))
        cp = c.lerp_cp((T0_in + T0_out) / 2, 0)

        # Compute the resulting efficiencies and power
        eta_s_tot = (pi_tot ** ((g - 1) / g) - 1) / (pi_tot ** ((g - 1) / (g * eta_p)) - 1)
        eta_s_stage = (pi_stage ** ((g - 1) / g) - 1) / (pi_stage ** ((g - 1) / (g * eta_p)) - 1)
        P = mdot * cp * T0_in * (1 / eta_s_tot * (pi_tot ** ((g - 1) / g) - 1))

        # Add the computed data for printing
        table.add_line(iter, 1e-6 * P, 100 * eta_s_tot, 100 * eta_s_stage, g)

        # Update gamma for next iteration
        g = cp / (cp - c.R_air)

    table.print()
