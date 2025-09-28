from matplotlib.pyplot import show
from argparse import ArgumentParser
import balloon_qnet.lowtran_piccia as lowtran
import sys

def horizontal(altitude, distance, wavelength, ihaze = 5):
    p = ArgumentParser(description = "Lowtran 7 interface")
    p.add_argument(
        "-z",
        "--obsalt",
        help = "altitude of both observers on horizontal path [km]",
        type = float,
        default = altitude,
    )
    p.add_argument(
        "-r",
        "--range_km",
        help = "range between observers on horizontal path [km]",
        type = float,
        default = distance,
    )
    p.add_argument(
        "-a",
        "--zenang",
        help="zenith angle [deg] can be single value or list of values",
        type = float,
        default = 0,
    )
    p.add_argument(
        "-i",
        "--ihaze",
        help="aerosol model",
        type = int,
        default = ihaze,
    )

    p.add_argument("-s", "--short", help = "shortest wavelength [nm]", type = float, default = 200)
    p.add_argument("-l", "--long", help = "longest wavelength [nm]", type = float, default = 30000)
    p.add_argument("-step", help = "wavelength step size [cm^-1]", type = float, default = 20)
    if "ipykernel" in sys.modules:
        P = p.parse_args(args=[])  # ignore Jupyter extra args
    else:
        P = p.parse_args()

    context = {
        "model": 5, # subarctic winter
        "itype": 1, # horizontal path
        "iemsct": 0, # tx mode
        "im": 0, # single scattering
        "ird1": 1, # card2C on
        "ihaze": P.ihaze, # urban aerosol
        "zmdl": P.obsalt, # altitude layer
        "h1": P.obsalt, # height of tx and rx
        "range_km": P.range_km, # distance between tx and rx
        "angle": P.zenang,
        "wlshort": P.short,
        "wllong": P.long,
        "wlstep": P.step,
    }

    TR = lowtran.lowtran(context).squeeze()

    return TR.sel(wavelength_nm = [wavelength], method = "nearest")["transmission"].values[0]

def slant(altitude_1, altitude_2, wavelength, zenith_angle, ihaze = 5):
    p = ArgumentParser(description = "Lowtran 7 interface")
    p.add_argument(
        "-z1",
        "--obsalt1",
        help = "altitude of one observer on slant path [km]",
        type = float,
        default = altitude_1,
    )
    p.add_argument(
        "-z2",
        "--obsalt2",
        help = "altitude of second observer on slant path [km]",
        type = float,
        default = altitude_2,
    )
    p.add_argument(
        "-a",
        "--zenang",
        help="zenith angle [deg] can be single value or list of values",
        type = float,
        default = zenith_angle,
    )
    p.add_argument(
        "-i",
        "--ihaze",
        help="aerosol model",
        type = int,
        default = ihaze,
    )
    p.add_argument("-s", "--short", help = "shortest wavelength [nm]", type = float, default = 200)
    p.add_argument("-l", "--long", help = "longest wavelength [nm]", type = float, default = 30000)
    p.add_argument("-step", help = "wavelength step size [cm^-1]", type = float, default = 20)
    if "ipykernel" in sys.modules:
        P = p.parse_args(args=[])  # ignore Jupyter extra args
    else:
        P = p.parse_args()

    context = {
        "model": 5, # subarctic winter
        "itype": 2, # slant path
        "iemsct": 0, # tx mode
        "im": 0, # single scattering
        "ird1": 1, # card2C on
        "ihaze": P.ihaze, # urban aerosol
        "h1": P.obsalt1, # height of observer 1
        "h2": P.obsalt2, # height of observer 2
        "angle": P.zenang,
        "wlshort": P.short,
        "wllong": P.long,
        "wlstep": P.step,
    }

    TR = lowtran.lowtran(context).squeeze()

    return TR.sel(wavelength_nm = [wavelength], method = "nearest")["transmission"].values[0]

