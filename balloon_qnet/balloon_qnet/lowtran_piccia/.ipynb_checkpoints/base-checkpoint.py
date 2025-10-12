from __future__ import annotations
from pathlib import Path
import importlib.util
import distutils.sysconfig
import logging
import xarray
import numpy as np
from typing import Any

from .cmake import build


def import_f2py_mod(name: str):

    mod_name = name + distutils.sysconfig.get_config_var("EXT_SUFFIX")  # type: ignore
    mod_file = Path(__file__).parent / mod_name
    if not mod_file.is_file():
        raise ModuleNotFoundError(mod_file)
    spec = importlib.util.spec_from_file_location(name, mod_file)
    if spec is None:
        raise ModuleNotFoundError(f"{name} not found in {mod_file}")
    mod = importlib.util.module_from_spec(spec)
    if mod is None:
        raise ImportError(f"could not import {name} from {mod_file}")
    spec.loader.exec_module(mod)  # type: ignore

    return mod


def nm2lt7(short_nm: float, long_nm: float, step_cminv: float = 20) -> tuple[float, float, float]:
    """
    converts wavelength in nm to cm^-1
    minimum meaningful step is 20, but 5 is minimum before crashing lowtran

    short: shortest wavelength e.g. 200 nm
    long: longest wavelength e.g. 30000 nm
    step: step size in cm^-1 e.g. 20

    output in cm^-1
    """
    short = 1e7 / short_nm
    long = 1e7 / long_nm

    N = int(np.ceil((short - long) / step_cminv)) + 1  # yes, ceil

    return short, long, N


def loopangle(context: dict[str, Any]) -> xarray.Dataset:
    """
    loop over "ANGLE" when context["angle"] is a vector
    """
    angles = np.atleast_1d(context["angle"])
    TR = xarray.Dataset(coords={"wavelength_nm": None, "angle_deg": angles})

    for a in angles:
        c = context.copy()
        c["angle"] = a
        TR = TR.merge(lowtran(c))

    return TR


def lowtran(context: dict[str, Any]) -> xarray.Dataset:
    """
    directly run Fortran code
    """
    # default parameters
    context.setdefault("time", None)
    context.setdefault("ird1", 0)
    context.setdefault("zmdl", 0)
    context.setdefault("p", 0)
    context.setdefault("t", 0)
    context.setdefault("wmol", [0] * 12)
    context.setdefault("h1", 0)
    context.setdefault("h2", 0)
    context.setdefault("iseasn", 0)
    context.setdefault("ivulc", 0)
    context.setdefault("icstl", 1)
    context.setdefault("icld", 0)
    context.setdefault("range_km", 0)
    

    # input check
    assert len(context["wmol"]) == 12, "see Lowtran user manual for 12 values of WMOL"
    assert np.isfinite(context["h1"]), "per Lowtran user manual Table 14, H1 must always be defined"
    # setup wavelength
    context.setdefault("wlstep", 20)
    if context["wlstep"] < 5:
        logging.critical("minimum resolution 5 cm^-1, specified resolution 20 cm^-1")

    wlshort, wllong, nwl = nm2lt7(context["wlshort"], context["wllong"], context["wlstep"])

    if not 0 < wlshort and wllong <= 50000:
        logging.critical("specified model range 0 <= wavelength [cm^-1] <= 50000")
        
    # invoke lowtran
    try:
        lowtran7 = import_f2py_mod("lowtran7")
    except ImportError:
        build(source_dir=Path(__file__).parent, build_dir=Path(__file__).parent / "build")
        lowtran7 = import_f2py_mod("lowtran7")

    Tx, V, Alam, trace, unif, suma, irrad, sumvv = lowtran7.lwtrn7(
        True,               # Enable Python interface
        nwl,                # wavelength step
        wllong,             # wavelength max
        wlshort,            # wavelength min
        context["wlstep"],  # wavelenght step in cm-1
        context["model"],   # card1 atmosphere model (0-7)
        context["itype"],   # card1 path type (1-3)
        context["iemsct"],  # card1 execution type (0-3)
        context["im"],      # card1 scattering off/on (0-1)
        context["ihaze"],   # card2 aerosol type (0-10)
        context["iseasn"],  # card2 seasonal aerosol (0-2)
        context["ivulc"],   # card2 aerosol profile and stratospheric aerosol (0-8)
        context["icstl"],   # card2 air mass character (1-10) >> HAZE 3 ONLY
        context["icld"],    # card2 cloud and rain models (0-20)
        context["ird1"],    # activate optional card2C off/on (0-1) for custom path
        context["zmdl"],    # card2C altitude layer boudnary [km]
        context["p"],       # card2C pressure layer boundary
        context["t"],       # card2C temperature layer boundary
        context["wmol"],    # card2C individual molecular species, see T11A p32
        context["h1"],      # card3 initial altitude [km]
        context["h2"],      # card3 final altitude [km] >> ITYPE 2 ONLY
        context["angle"],   # card3 initial zenith angle from h1 [deg]
        context["range_km"],# card3 path length [km]
    )
    
    dims = ("time", "wavelength_nm", "angle_deg")
    TR = xarray.Dataset(
        {
            "transmission": (dims, Tx[:, 9][None, :, None]),
            "radiance": (dims, sumvv[None, :, None]),
            "irradiance": (dims, irrad[:, 0][None, :, None]),
            "pathscatter": (dims, irrad[:, 2][None, :, None]),
        },
        coords={
            "time": [context["time"]],
            "wavelength_nm": Alam * 1e3,
            "angle_deg": [context["angle"]],
        },
    )

    return TR
