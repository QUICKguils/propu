import numpy as np

from propu import _SRC_PATH

_DATA_PATH = _SRC_PATH / "project" / "res"


def main() -> None:
    data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_11x10_geom.txt"), skiprows=1)
    span = 11 / 2
    stations = data[:, 0] * span
    pitches = np.deg2rad(data[:, 2])
    pitches_inch = 2 * np.pi * stations * np.tan(pitches)
    print(pitches_inch)
    print(np.mean(pitches_inch))


if __name__ == "__main__":
    main()
