import math
import numpy as np
import netsquid as ns
from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import DownlinkChannel, CachedChannel

# Parameters (you can adjust these)
wavelength = 1550e-9
ground_station_alt = 0 # 0.020
obs_ratio_ground = 0 # 0.3
Cn0 = 9.6e-14
u_rms = 10
pointing_error = 1e-6
Qonnector_meas_succ = 0.85
tracking_efficiency = 1 # 0.8
rx_aperture_ground = 0.3
W0 = 0.1
simtime = 100000   # shorter for testing

#Parameters and function to calculate the secret key rate
ratesources = 80e6
sourceeff = 0.01
QBER = 0.04

def theoretical_eff(distance, h_balloons, n):
    ## Assuming flat earth approximation
    zenith_angle = np.rad2deg(np.pi/2 - np.arcsin(h_balloons/distance))
    
    Tatm = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    ch = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n, Cn0, u_rms,
                         wavelength, ground_station_alt, h_balloons, zenith_angle,
                         pointing_error=pointing_error, tracking_efficiency=tracking_efficiency, Tatm=Tatm)
    eta = np.arange(1e-7, 1, 0.001)
    return ch._compute_mean_channel_efficiency(eta, distance, detector_efficiency=Qonnector_meas_succ)





def simulated_eff_requires_fix(distance, h_balloons, n):
    ## Assuming flat earth approximation
    zenith_angle = np.rad2deg(np.pi/2 - np.arcsin(h_balloons/distance))
    
    Tatm = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    ch = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n, Cn0, u_rms,
                         wavelength, ground_station_alt, distance, zenith_angle,
                         pointing_error=pointing_error, tracking_efficiency=tracking_efficiency, Tatm=Tatm)
    
    a = ch._compute_loss_probability(distance, math.ceil(simtime/Qonnector_init_time))
    net = QEurope("Europe")
    net.Add_Qonnector("Alice")
    net.Add_Qonnector("Bob")
    downlink = CachedChannel(a)
    net.connect_qonnectors("Alice", "Bob", distance=distance, loss_model=downlink)
    
    alice = net.network.get_node("Alice")
    bob = net.network.get_node("Bob")
    send = SendBB84(alice, Qonnector_init_succ, Qonnector_init_flip, bob)
    receive = ReceiveProtocol(bob, Qonnector_meas_succ, Qonnector_meas_flip, True, alice)
    send.start(); receive.start()
    ns.sim_run(duration=simtime)
    
    if len(bob.QlientKeys[alice.name]) == 0: return 0.0
    return len(alice.QlientKeys[bob.name]) / len(bob.QlientKeys[alice.name])








def simulated_eff(distance, h_balloons, n, simtime=100000):
    ## Assuming flat earth approximation
    zenith_angle = np.rad2deg(np.pi/2 - np.arcsin(h_balloons/distance))
    
    Tatm = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    ch = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n, Cn0, u_rms,
                         wavelength, ground_station_alt, distance, zenith_angle,
                         pointing_error=pointing_error, tracking_efficiency=tracking_efficiency, Tatm=Tatm)
    
    # Probability of loss per qubit
    loss_probs = ch._compute_loss_probability(distance, simtime)
    
    # Monte Carlo simulate
    sent = simtime
    received = np.random.binomial(sent, 1-loss_probs.mean())  # approximate survival
    return received / sent

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
        Secret key rate in kbit/s.
    """
    if efficiency <= 0:
        return 0.0
    rate = ratesources * sourceeff * efficiency
    skr = rate * (1 - 2 * h(QBER)) / 1000
    return skr


# distance = 25  # km (LoS distance)
# n = 5          # AO correction order

# theo = theoretical_eff(distance, 15, n)
# simu = simulated_eff(distance, 15, n)

# print("Theoretical mean efficiency:", theo)
# print("Simulated mean efficiency:", simu)