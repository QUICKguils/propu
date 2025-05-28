from numpy import *
import pandas
import scipy

# ------------------------------------------------------------------------------
# read all available characteristic data
# ------------------------------------------------------------------------------
speeds = [90, 95, 100]
perfo = {}
for s in speeds:
    fname = "speedLine_%sNn.csv" % s
    perfo[s] = pandas.read_csv(
        fname, dtype=float64, header=None, names=["mFlow", "pRatio", "efficiency"]
    )


# ------------------------------------------------------------------------------
# interpolate (pressure ratio, polytropic efficiency) ifo mass flow on characteristic
# - perfo : characteristic
# - mFlow : fraction of nominal corrected mass flow rate
# ------------------------------------------------------------------------------
def performance1(perfoNn, mFlow):
    p = interp(mFlow, perfoNn["mFlow"], perfoNn["pRatio"])
    e = interp(mFlow, perfoNn["mFlow"], perfoNn["efficiency"])
    return (p, e)


# ------------------------------------------------------------------------------
# find (pressure ratio, polytropic efficiency) ifo speed and mass flow
# - speed : percentage (integer) of nominal corrected rotation speed
# - mFlow : fraction of nominal corrected mass flow rate
# ------------------------------------------------------------------------------
def performance(speed, mFlow):
    return performance1(perfo[speed], mFlow)
