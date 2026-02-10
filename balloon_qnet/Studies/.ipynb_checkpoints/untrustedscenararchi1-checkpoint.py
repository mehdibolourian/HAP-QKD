from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import DownlinkChannel, CachedChannel, RE

"""This script estimates and print the Entanglement-based QKD secret key rate between the two Qlients of the Italian network,
 as described in scenario of Fig20(a) of the paper """

# Parameters
wavelength = 1550e-9
EPR_succ = 1
p_transmit=1
ground_station_alt = 0.020 #Altitude of the receiving telescopes
W0 = 0.1 #Initial Beam Waist 
obs_ratio_ground = 0.3 #Obscuration ratio of the receiving telescope same
n_max_ground = 6 #Maximum radial index of correction of AO system 
Cn0 = 9.6*10**(-14) #Reference index of refraction structure constant at ground level
u_rms = 10 #Wind speed
pointing_error = 1e-6 #Pointing error variance
Qonnector_meas_succ = 0.85 #Detector efficiency at the receiver
tracking_efficiency = 0.8 #
h_balloons = 35 #Altitude of the balloon (km)
rx_aperture = 0.4 #Aperture of the receiving telescope (cm)
dist_cities = 189 #distance between the cities (km)
distAlice = 36 #distance between Alice and her Qonnector (km)
distBob = 50 #distance between Bob and her Qonnector (km)

#Parameter to estimate the secret key rate
ratesources = 80e6
sourceeff=0.01
QBER = 0.04

def h(p):
    """Binary entropy function"""
    return -p*np.log2(p)-(1-p)*np.log2(1-p)


simtime = 2000000 #Simulation time

def GetRate(h_balloons,dist_cities):
    ns.sim_reset()
    alpha = dist_cities/(2*RE)
    RG = RE + ground_station_alt
    RS = RE + h_balloons
    chan_length = np.sqrt(RG**2 + RS**2 - 2*RG*RS*np.cos(alpha))
    print("channel length: "+ str(chan_length))
    zenith_angle = np.rad2deg(np.arccos((RS**2 - RG**2 - chan_length**2 )/(2*chan_length*RG)))
    print("zenith angle: " + str(zenith_angle))
    transmittance_vert = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    # Initialize network
    net = QEurope("Europe")

    # Create two quantum Cities
    net.Add_Qonnector("QonnectorPadova")
    net.Add_Qlient("QlientAlice",distAlice,"QonnectorPadova")

    net.Add_Qonnector("QonnectorFlorence")
    net.Add_Qlient("QlientBob",distBob,"QonnectorFlorence")
    
    # Create drone 1
    net.Add_Qonnector("QonnectorDrone")
    
    # Create channels
    downlink_channel1 = DownlinkChannel(W0, rx_aperture, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, 
                                                          h_balloons, zenith_angle, pointing_error = pointing_error, tracking_efficiency = tracking_efficiency,
                                                            Tatm = transmittance_vert)
    downlink_channel2 = DownlinkChannel(W0, rx_aperture, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, 
                                                          h_balloons, zenith_angle, pointing_error = pointing_error, tracking_efficiency = tracking_efficiency,
                                                            Tatm = transmittance_vert)
    
    a = downlink_channel1._compute_loss_probability(chan_length,math.ceil(simtime/Qonnector_init_time ))
    b = downlink_channel2._compute_loss_probability(chan_length,math.ceil(simtime/Qonnector_init_time ))
    downlink1 = CachedChannel(a)
    downlink2 = CachedChannel(b)
    # Connect nodes
    net.connect_qonnectors("QonnectorPadova", "QonnectorDrone", distance = chan_length, loss_model = downlink1)
    net.connect_qonnectors("QonnectorFlorence", "QonnectorDrone", distance = chan_length, loss_model = downlink2)

    # Get node instances
    city1 = net.network.get_node("QonnectorPadova")
    city2 = net.network.get_node("QonnectorFlorence")
    balloon = net.network.get_node("QonnectorDrone")
    Alice = net.network.get_node("QlientAlice")
    Bob= net.network.get_node("QlientBob")

    # Entanglement based QDK from city 1 to drone 1
    send = SendEPR(city1, city2, EPR_succ, balloon)
    send.start()
    
    transmit1 = TransmitProtocol(balloon, Alice, p_transmit, city1)
    transmit1.start()

    transmit2 = TransmitProtocol(balloon, Bob, p_transmit, city2)
    transmit2.start()

    receive1 = ReceiveProtocol(city1, Qonnector_meas_succ, Qonnector_meas_flip, True, Alice)
    receive1.start()
    receive2 = ReceiveProtocol(city2, Qonnector_meas_succ, Qonnector_meas_flip, True, Bob)
    receive2.start()

    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display results 
    L1 = Sifting(Alice.keylist, Bob.keylist) 
    sent = len(balloon.QlientKeys[city1.name])
    rec = len(L1)
    eff = rec/sent
    rate = ratesources*sourceeff*eff
    skr = rate*(1-2*h(QBER))

    print("Channel length: "+str(chan_length)+"km")
    print("Zenith angle: " + str(zenith_angle) + ' degree')
    print("EPR pairs sent from balloon: " + str(sent))
    print("Number of qubits received in both city: " +str(rec) )
    print("Channel efficiency : "+str(eff) + " bits per channel use")
    print("Secret key rate: "+ str(skr) + " bits per second")

    return eff,rate, skr


def AverageRates(n):
    toteff = []
    totrate = []
    totskr = []
    for i in range(n):
        print("test n"+str(i))
        eff,rate,skr = GetRate(h_balloons,dist_cities)
        toteff.append(eff)
        totrate.append(rate)
        totskr.append(skr)
    return toteff, totrate, totskr

toteff, totrate, totskr = AverageRates(5)  
eff = np.mean(toteff)
effstd = np.std(toteff)
rate = np.mean(totrate)
ratestd = np.std(totrate)
skr = np.mean(totskr)
skrstd = np.std(totskr)

print("E91 rates  \n")
print("Alice and Bob: \n")
print("Rate: " + str(rate) + "+/-" + str(ratestd) + "\n")
print("SKR: " + str(skr) + "+/-" + str(skrstd) + "\n")
print("Chann efficiency: " + str(eff) + "+/-" + str(effstd) + "\n")