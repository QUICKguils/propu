from propu.util import isatmosphere
from propu.util.constant import ur

def main():
    print("Solving exo 2.2: Sonaca 200")

    # Statement data
    n_blades = 3
    diameter = (1750*ur.mm).to(ur.m)
    rho, p, T, a = isatmosphere(0*ur.m)
    thrust = 1700*ur.N
    v_inf = (150*ur("km/hr")).to("m/s")

    print(f"{v_inf=:~P}p")
    print(f"{p=:~P}")
