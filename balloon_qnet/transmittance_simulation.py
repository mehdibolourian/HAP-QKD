import math
import numpy as np
import netsquid as ns
import balloon_qnet.transmittance as transmittance
# from balloon_qnet.QEuropeFunctions import *
from balloon_qnet.free_space_losses import (
    UplinkChannel,
    DownlinkChannel,
    compute_channel_length,
    CachedChannel
)

# Parameters (you can adjust these)
wavelength = 1550e-9
ground_station_alt = 0 # 0.020
obs_ratio_ground = 0.3
Cn0 = 9.6e-14
u_rms = 10
pointing_error = 1e-6
Qonnector_meas_succ = 0.85
tracking_efficiency = 0.8
rx_aperture_ground = 0.6 # 0.3
W0 = 0.1
simtime = 500   # shorter for testing

#Parameters and function to calculate the secret key rate
ratesources = 80e6
sourceeff = 0.01
QBER = 0.04

params = {
    "wavelength": 1550e-9,
    "W0": 0.1,
    "rx_aperture_down": 0.6,
    "rx_aperture_up": 0.6,
    "obs_ratio": 0.3,
    "Cn0": 9.6e-14,
    "u_rms": 10,
    "pointing_error": 1e-6,
    "tracking_efficiency": 0.8,
    "detector_eff": 0.85,
    "init_time": 1  # example value
}

# ---------------------------
# Helper: compute zenith angle from LoS distance and balloon height
# ---------------------------
def compute_zenith_angle(gs_alt, balloon_alt, distance):
    """
    Compute zenith angle (radians) given:
      gs_alt: ground station altitude (km)
      balloon_alt: balloon altitude (km)
      distance: line-of-sight distance between ground and balloon (km)
    """
    vertical = balloon_alt - gs_alt
    # Clamp to avoid invalid acos values
    cos_theta = max(-1.0, min(1.0, vertical / distance))
    return np.degrees(math.acos(cos_theta))


# ---------------------------
# Theoretical mean efficiency
# ---------------------------
def channel_theory(direction, gs_alt, balloon_alt, distance, n_correction):
    """
    direction: "uplink" or "downlink"
    gs_alt: ground station altitude (km)
    balloon_alt: balloon altitude (km)
    distance: line-of-sight distance (km)
    n_correction: AO system correction index
    """
    zenith_angle = compute_zenith_angle(gs_alt, balloon_alt, distance)

    Tatm = transmittance.slant(gs_alt, balloon_alt, params["wavelength"] * 1e9, zenith_angle)

    print(f"Tatm: {Tatm}")

    if direction == "uplink":
        channel = UplinkChannel(
            params["W0"], params["rx_aperture_up"], params["obs_ratio"],
            n_correction, params["Cn0"], params["u_rms"],
            params["wavelength"], gs_alt, balloon_alt, zenith_angle,
            pointing_error=params["pointing_error"],
            tracking_efficiency=params["tracking_efficiency"],
            Tatm=Tatm
        )
    elif direction == "downlink":
        channel = DownlinkChannel(
            params["W0"], params["rx_aperture_down"], params["obs_ratio"],
            n_correction, params["Cn0"], params["u_rms"],
            params["wavelength"], gs_alt, balloon_alt, zenith_angle,
            pointing_error=params["pointing_error"],
            tracking_efficiency=params["tracking_efficiency"],
            Tatm=Tatm
        )
    else:
        raise ValueError("direction must be 'uplink' or 'downlink'")

    eta = np.arange(1e-7, 1, 0.001)
    mean = channel._compute_mean_channel_efficiency(
        eta, distance, detector_efficiency=params["detector_eff"]
    )
    return mean


# ---------------------------
# Simulated efficiency
# ---------------------------
def channel_simulation(direction, gs_alt, balloon_alt, distance, n_correction, simtime=simtime):
    """
    direction: "uplink" or "downlink"
    gs_alt, balloon_alt, distance: geometry (km)
    params: dict of channel parameters
    n_correction: AO correction index
    simtime: number of channel uses to simulate
    """
    zenith_angle = zenith_angle = compute_zenith_angle(0, balloon_alt, distance)

    Tatm = transmittance.slant(gs_alt, balloon_alt, params["wavelength"] * 1e9, zenith_angle)

    if direction == "uplink":
        channel = UplinkChannel(
            params["W0"], params["rx_aperture_up"], params["obs_ratio"],
            n_correction, params["Cn0"], params["u_rms"],
            params["wavelength"], gs_alt, balloon_alt, zenith_angle,
            pointing_error=params["pointing_error"],
            tracking_efficiency=params["tracking_efficiency"],
            Tatm=Tatm
        )
    elif direction == "downlink":
        channel = DownlinkChannel(
            params["W0"], params["rx_aperture_down"], params["obs_ratio"],
            n_correction, params["Cn0"], params["u_rms"],
            params["wavelength"], gs_alt, balloon_alt, zenith_angle,
            pointing_error=params["pointing_error"],
            tracking_efficiency=params["tracking_efficiency"],
            Tatm=Tatm
        )
    else:
        raise ValueError("direction must be 'uplink' or 'downlink'")

    # Simulate by sampling loss probabilities
    a = channel._compute_loss_probability(distance, math.ceil(simtime / params["init_time"]))

    # Channel efficiency = fraction of photons surviving
    efficiency = 1.0 - np.mean(a)  # or another stat depending on how a is defined
    
    # Monte Carlo simulate
    # sent = simtime / params["init_time"]
    # received = np.random.binomial(sent, 1 - a.mean())  # approximate survival
    # efficiency = received / sent

    # print("a type:", type(a), "len:", len(a))
    # a = np.asarray(a)
    # print("a dtype:", a.dtype, "min,max,mean:", a.min(), a.max(), a.mean())
    # print("first 20:", a[:20])
    return efficiency





# def theoretical_eff(distance, h_balloons, n):
#     ## Assuming flat earth approximation
#     zenith_angle = np.rad2deg(np.pi/2 - np.arcsin(h_balloons/distance))
    
#     Tatm = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
#     ch = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n, Cn0, u_rms,
#                          wavelength, ground_station_alt, h_balloons, zenith_angle,
#                          pointing_error=pointing_error, tracking_efficiency=tracking_efficiency, Tatm=Tatm)
#     eta = np.arange(1e-7, 1, 0.001)
#     return ch._compute_mean_channel_efficiency(eta, distance, detector_efficiency=Qonnector_meas_succ)


# def simulated_eff(distance, h_balloons, n, simtime=simtime):
#     ## Assuming flat earth approximation
#     zenith_angle = zenith_angle = compute_zenith_angle(0, h_balloons, distance) #np.rad2deg(np.pi/2 - np.arcsin(h_balloons/distance))
    
#     Tatm = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    
#     ch = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n, Cn0, u_rms,
#                          wavelength, ground_station_alt, h_balloons, zenith_angle,
#                          pointing_error=pointing_error, tracking_efficiency=tracking_efficiency, Tatm=Tatm)
    
#     # Probability of loss per qubit
#     loss_probs = ch._compute_loss_probability(distance, simtime)
    
#     # Monte Carlo simulate
#     sent = simtime
#     received = np.random.binomial(sent, 1-loss_probs.mean())  # approximate survival
#     return received / sent

def h(p):
    """Binary entropy function"""
    return -p*np.log2(p)-(1-p)*np.log2(1-p)

def compute_skr(efficiency, ratesources=80e6, sourceeff=0.01, QBER=0.04):
    """
    Compute BB84 secret key rate (kbit/s) given channel efficiency.

    Parameters
    ----------
    efficiency : float
        Channel transmission efficiency (0–1), from theory or simulation.
    ratesources : float
        Pulse repetition rate (Hz).
    sourceeff : float
        Source efficiency.
    QBER : float
        Quantum Bit Error Rate (0–1).

    Returns
    -------
    skr : float
        Secret key rate in bit/s.
    """
    if efficiency <= 0:
        return 0.0
    rate = ratesources * sourceeff * efficiency
    skr = rate * (1 - 2 * h(QBER))
    return skr