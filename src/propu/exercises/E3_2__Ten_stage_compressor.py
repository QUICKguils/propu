import propu.constant as c


class Table:
    """Iteration table printer."""

    hline = "-" * 53
    variables = f"|{'Iter':^6}|{'Power':^11}|{'eta_s_tot':^11}|{'eta_s_stage':^11}|{'gamma':^8}|"
    units = f"|{'':^6}|{'(MW)':^11}|{'(%)':^11}|{'(%)':^11}|{'':^8}|"

    header = (hline, variables, units, hline)
    footer = hline

    def __init__(self):
        self.lines = []

    def add_line(self, iter, P, eta_s, eta_p, g):
        self.lines.append(
            f"|{iter:^6}|{1e-6 * P:^11.4f}|{100 * eta_s:^11.2f}|{100 * eta_p:^11.2f}|{g:^8.4f}|"
        )

    def print(self):
        print(*Table.header, sep="\n")
        print(*self.lines, sep="\n")
        print(Table.footer)


def main():
    print("Solving exo 3.2: Ten-stage compressor")

    # Statement data
    mdot = 164 * c.uconv("lb/s", "kg/s")
    pi_stage = 1.3
    n_stage = 10
    eta_p = 0.9
    p0_in = 57_452  # [Pa]
    T0_in = 266.42  # [K]

    # Instantiate a table printer
    table = Table()

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

        # Update gamma for next iteration
        g = cp / (cp - c.R_air)

        # Add the computed data for printing
        table.add_line(iter, P, eta_s_tot, eta_s_stage, g)

    table.print()
