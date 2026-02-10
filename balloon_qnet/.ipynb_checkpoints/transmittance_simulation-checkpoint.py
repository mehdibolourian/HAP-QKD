import math
import numpy as np
import netsquid as ns

import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import (
    UplinkChannel,
    DownlinkChannel,
    CachedChannel
)
from balloon_qnet.QEuropeFunctions import *

# ============================================================
# Physical constants
# ============================================================

RE = 6371.0  # Earth radius [km]

# ============================================================
# Channel & system parameters (paper-consistent)
# ============================================================


# "W0": 0.15,
#     "rx_aperture_down": 1.0,
#     "rx_aperture_up": 0.50,
#     "pointing_error": 2e-6,
#     "tracking_efficiency": 0.9,
params = {
    # Optical
    "wavelength": 1550e-9,     # [m]
    "W0": 0.15,                # Initial beam waist [m] (15 cm)
    
    # Apertures
    "rx_aperture_down": 0.75,  # Ground station [m]
    "rx_aperture_up": 0.40,    # HAP [m]
    "obs_ratio": 0.3,

    # Atmosphere & turbulence
    "Cn0": 1e-13,            # [m^(-2/3)] "Cn0": 9.6e-14,            # [m^(-2/3)]
    "u_rms": 10,           # [m/s]

    # Pointing & tracking
    "pointing_error": 1e-7,    # 1 μrad (paper)
    "tracking_efficiency": 0.85,

    # Detection
    "detector_eff": 0.85,

    # Netsquid
    "init_time": 1
}

# ============================================================
# BB84 / SKR parameters
# ============================================================

ratesources = 80e6
sourceeff = 0.01
QBER = 0.04
simtime = 50_000

# ============================================================
# Geometry helpers (spherical Earth)
# ============================================================

def slant_range(gs_alt, hap_alt, zenith_angle_deg):
    """
    Compute slant range using spherical Earth geometry.

    gs_alt, hap_alt : km
    zenith_angle_deg: degrees (at GS)
    """
    theta = np.radians(zenith_angle_deg)
    RG = RE + gs_alt
    RH = RE + hap_alt

    return np.sqrt(
        RH**2 + RG**2 * (np.cos(theta)**2 - 1) - RG * np.cos(theta)
    )


def zenith_angle_from_distance(gs_alt, hap_alt, distance):
    """
    Inverse: compute zenith angle from slant distance.
    """
    RG = RE + gs_alt
    RH = RE + hap_alt
    L = distance

    cos_theta = -(RG**2 + L**2 - RH**2) / (2 * RG * L)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))

# ============================================================
# Theoretical mean channel efficiency
# ============================================================

def channel_theory(direction, gs_alt, hap_alt, distance, n_correction, params_in):
    """
    Compute theoretical mean channel efficiency.
    """

    zenith = zenith_angle_from_distance(gs_alt, hap_alt, distance)

    Tatm = transmittance.slant(
        gs_alt,
        hap_alt,
        params_in["wavelength"] * 1e9,
        zenith
    )

    if direction == "uplink":
        channel = UplinkChannel(
            params_in["W0"],
            params_in["rx_aperture_up"],
            params_in["obs_ratio"],
            n_correction,
            params_in["Cn0"],
            params_in["u_rms"],
            params_in["wavelength"],
            gs_alt,
            hap_alt,
            zenith,
            pointing_error=params_in["pointing_error"],
            tracking_efficiency=params_in["tracking_efficiency"],
            Tatm=Tatm
        )

    elif direction == "downlink":
        channel = DownlinkChannel(
            params_in["W0"],
            params_in["rx_aperture_down"],
            params_in["obs_ratio"],
            n_correction,
            params_in["Cn0"],
            params_in["u_rms"],
            params_in["wavelength"],
            gs_alt,
            hap_alt,
            zenith,
            pointing_error=params_in["pointing_error"],
            tracking_efficiency=params_in["tracking_efficiency"],
            Tatm=Tatm
        )
    else:
        raise ValueError("direction must be 'uplink' or 'downlink'")

    eta = np.arange(1e-7, 1.0, 1e-3)
    mean_eta = channel._compute_mean_channel_efficiency(
        eta,
        distance,
        detector_efficiency=params_in["detector_eff"]
    )

    return mean_eta

# ============================================================
# Monte Carlo channel simulation (NetSquid)
# ============================================================

def channel_simulation(direction, gs_alt, hap_alt, distance, n_correction):
    """
    Simulate BB84 over free-space HAP–GS channel.
    """

    zenith = zenith_angle_from_distance(gs_alt, hap_alt, distance)

    Tatm = transmittance.slant(
        gs_alt,
        hap_alt,
        params["wavelength"] * 1e9,
        zenith
    )

    # Initialize network
    net = QEurope("HAPNet")
    net.Add_Qonnector("GS")
    net.Add_Qonnector("HAP")

    if direction == "uplink":
        channel = UplinkChannel(
            params["W0"],
            params["rx_aperture_up"],
            params["obs_ratio"],
            n_correction,
            params["Cn0"],
            params["u_rms"],
            params["wavelength"],
            gs_alt,
            hap_alt,
            zenith,
            pointing_error=params["pointing_error"],
            tracking_efficiency=params["tracking_efficiency"],
            Tatm=Tatm
        )
    elif direction == "downlink":
        channel = DownlinkChannel(
            params["W0"],
            params["rx_aperture_down"],
            params["obs_ratio"],
            n_correction,
            params["Cn0"],
            params["u_rms"],
            params["wavelength"],
            gs_alt,
            hap_alt,
            zenith,
            pointing_error=params["pointing_error"],
            tracking_efficiency=params["tracking_efficiency"],
            Tatm=Tatm
        )
    else:
        raise ValueError("direction must be 'uplink' or 'downlink'")

    # Sample channel loss
    loss_samples = channel._compute_loss_probability(
        distance,
        math.ceil(simtime / params["init_time"])
    )

    loss_model = CachedChannel(loss_samples)

    # Connect GS <-> HAP with correct distance
    net.connect_qonnectors(
        "GS",
        "HAP",
        distance=distance,
        loss_model=loss_model
    )

    gs = net.network.get_node("GS")
    hap = net.network.get_node("HAP")

    send = SendBB84(hap, Qonnector_init_succ, Qonnector_init_flip, gs)
    recv = ReceiveProtocol(gs, params["detector_eff"], Qonnector_meas_flip, True, hap)

    send.start()
    recv.start()

    ns.sim_run(simtime)

    sifted = Sifting(gs.QlientKeys[hap.name], hap.QlientKeys[gs.name])
    efficiency = len(hap.QlientKeys[gs.name]) / max(1, len(gs.QlientKeys[hap.name]))

    return efficiency

# ============================================================
# Secret key rate
# ============================================================

def h(p):
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


def compute_skr(efficiency):
    if efficiency <= 0:
        return 0.0
    raw_rate = ratesources * sourceeff * efficiency
    
    return raw_rate * (1 - 2 * h(QBER))
