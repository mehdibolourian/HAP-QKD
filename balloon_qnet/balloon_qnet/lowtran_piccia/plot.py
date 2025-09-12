from __future__ import annotations
import numpy as np
import xarray
from matplotlib.pyplot import figure
from typing import Any

#
h = 6.62607004e-34 # Planck's constant
c = 299792458 # speed of light
UNITS = r"ster$^{-1}$ cm$^{-2}$ $\mu$m$^{-1}$]"
plotNp = False


def scatter(irrad: xarray.Dataset, context: dict[str, Any], log: bool = False) -> None:

    fg = figure()
    axs = fg.subplots(2, 1, sharex=True)

    transtxt = "Transmittance"

    ax = axs[0]
    ax.plot(irrad.wavelength_nm, irrad["transmission"].squeeze())
    ax.set_title(transtxt)
    ax.set_ylabel("Transmission (unitless)")
    ax.grid(True)
    ax.legend(irrad.angle_deg.values)

    ax = axs[1]
    if plotNp:
        Np = (irrad["pathscatter"] * 10000) * (irrad.wavelength_nm * 1e9) / (h * c)
        ax.plot(irrad.wavelength_nm, Np)
        ax.set_ylabel("Photons [s$^{-1}$ " + UNITS)
    else:
        ax.plot(irrad.wavelength_nm, irrad["pathscatter"].squeeze())
        ax.set_ylabel("Radiance [W " + UNITS)

    ax.set_xlabel("wavelength [nm]")
    ax.set_title("Single-scatter Path Radiance")
    ax.invert_xaxis()
    ax.autoscale(True, axis="x", tight=True)
    ax.grid(True)

    if log:
        ax.set_yscale("log")
    #        ax.set_ylim(1e-8,1)

    try:
        fg.suptitle(f'Obs. to Space: zenith angle: {context["angle"]} deg., ')
        # {datetime.utcfromtimestamp(irrad.time.item()/1e9)}
    except (AttributeError, TypeError):
        pass


def radiance(irrad: xarray.Dataset, context: dict[str, Any], log: bool = False) -> None:
    fg = figure()
    axs = fg.subplots(2, 1, sharex=True)

    transtxt = "Transmittance Observer to Space"

    ax = axs[0]
    ax.plot(irrad.wavelength_nm, irrad["transmission"].squeeze())
    ax.set_title(transtxt)
    ax.set_ylabel("Transmission (unitless)")
    ax.grid(True)

    ax = axs[1]
    if plotNp:
        Np = (irrad["radiance"] * 10000) * (irrad.wavelength_nm * 1e9) / (h * c)
        ax.plot(irrad.wavelength_nm, Np)
        ax.set_ylabel("Photons [s$^{-1}$ " + UNITS)
    else:
        ax.plot(irrad.wavelength_nm, irrad["radiance"].squeeze())
        ax.set_ylabel("Radiance [W " + UNITS)

    ax.set_xlabel("wavelength [nm]")
    ax.set_title("Atmospheric Radiance")
    ax.invert_xaxis()
    ax.autoscale(True, axis="x", tight=True)
    ax.grid(True)

    if log:
        ax.set_yscale("log")
        ax.set_ylim(1e-8, 1)

    try:
        fg.suptitle(f'Obs. zenith angle: {context["angle"]} deg., ')
        # {datetime.utcfromtimestamp(irrad.time.item()/1e9)}
    except (AttributeError, TypeError):
        pass


def transmission(T: xarray.Dataset, context: dict[str, Any], log: bool = False) -> None:
    ax = figure().gca()

    h = ax.plot(T.wavelength_nm, T["transmission"].squeeze())

    ax.set_xlabel("wavelength [nm]")
    ax.set_ylabel("transmission (unitless)")
    ax.set_title(f'Transmittance Ground-Space: Obs. zenith angle {context["angle"]} deg.')
    # ax.legend(loc='best')
    ax.grid(True)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 1)
    else:
        ax.set_ylim(0, 1)
    ax.invert_xaxis()
    ax.autoscale(True, axis="x", tight=True)
    ax.legend(h, T.angle_deg.values)


def irradiance(irrad: xarray.Dataset, context: dict[str, Any], log: bool = False) -> None:
    fg = figure()
    axs = fg.subplots(2, 1, sharex=True)

    #    if context['isourc'] == 0:
    stxt = "Sun's"
    #    elif context['isourc'] == 1:
    #        stxt = "Moon's"
    #    else:
    #        raise ValueError(f'ISOURC={context["isourc"]} not defined case')

    stxt += f' zenith angle {irrad.angle_deg.values} deg., Obs. height {context["h1"]} km. '
    try:
        stxt += np.datetime_as_string(irrad.time)[:-10]
    except (AttributeError, TypeError):
        pass

    fg.suptitle(stxt)

    if context["iemsct"] == 3:
        key = "irradiance"
        transtxt = "Transmittance Observer to Space"
    elif context["iemsct"] == 1:
        key = "radiance"
        transtxt = "Transmittance Observer to Observer"

    # irrad.['transmission'].plot()

    ax = axs[0]
    h = ax.plot(irrad.wavelength_nm, irrad["transmission"].squeeze())
    ax.set_title(transtxt)
    ax.set_ylabel("Transmission (unitless)")
    ax.grid(True)
    try:
        ax.legend(h, irrad.angle_deg.values)
    except AttributeError:
        pass

    ax = axs[1]
    ax.plot(irrad.wavelength_nm, irrad[key].squeeze())
    ax.set_xlabel("wavelength [nm]")
    ax.invert_xaxis()
    ax.grid(True)

    if context["iemsct"] == 3:
        ttxt = "Irradiance "
        ax.set_ylabel("Solar Irradiance [W " + UNITS)
        ax.set_title(ttxt)
    elif context["iemsct"] == 1:
        ttxt = "Single-scatter Radiance "
        ax.set_ylabel("Radiance [W " + UNITS)
        ax.set_title(ttxt)

    if log:
        ax.set_yscale("log")
        ax.set_ylim(1e-8, 1)

    ax.autoscale(True, axis="x", tight=True)


def horiz(trans: xarray.Dataset, context: dict[str, Any], log: bool = False) -> None:

    ttxt = f'Transmittance Horizontal \n {context["range_km"]} km path @ {context["h1"]} km altitude\n'

    if context["model"] == 0:
        ttxt += f'User defined atmosphere: pressure: {context["p"]} mbar, temperature {context["t"]} K'
    elif context["model"] == 5:
        ttxt += "Subarctic winter atmosphere"

    ax = figure().gca()

    ax.plot(trans.wavelength_nm, trans["transmission"].squeeze())

    ax.set_xlabel("wavelength [nm]")
    ax.set_ylabel("transmission (unitless)")
    ax.set_title(ttxt)
    # ax.legend(loc='best')
    ax.grid(True)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 1)
    else:
        ax.set_ylim(0, 1)
    ax.invert_xaxis()
    ax.autoscale(True, axis="x", tight=True)
