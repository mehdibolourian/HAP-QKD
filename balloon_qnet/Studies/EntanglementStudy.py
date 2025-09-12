from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import DownlinkChannel, CachedChannel, RE
import multiprocessing as mlp
import os
import functools as fnct

'''This script calculates the number of successful EPR pairs generation per second between two cities as a function of the distance between them,
    using a balloon and free space link or using a ground-station and fiber links.
   It creates 2 files, EPRFreespace.txt and EPRFiber.txt, containing each the EPR transmission rates for respectively the balloon link and the fiber link
'''

# Parameters
wavelength = 1550e-9
EPR_succ = 1
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
p_transmit=1

dist_cities = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200] #distance between the cities (km)

#Parameters and function to calculate the secret key rate
ratesources = 80e6 
sourceeff=0.01
QBER = 0.04

def h(p):
    """Binary entropy function"""
    return -p*np.log2(p)-(1-p)*np.log2(1-p)


simtime = 10000000 #Simulation time

def GetRatefree(dist_cities):
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

    net.Add_Qonnector("QonnectorFlorence")
    
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

    # Entanglement based QDK from city 1 to drone 1
    send = SendEPR(city1, city2, EPR_succ, balloon)
    send.start()

    receive1 = ReceiveProtocol(balloon, Qonnector_meas_succ, Qonnector_meas_flip, True, city1)
    receive1.start()
    receive2 = ReceiveProtocol(balloon, Qonnector_meas_succ, Qonnector_meas_flip, True, city2)
    receive2.start()

    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display results 
    L1 = Sifting(city1.QlientKeys[balloon.name], city2.QlientKeys[balloon.name])
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

def GetRatefiber(dist_cities):
    ns.sim_reset()

    # Initialize network
    net = QEurope("Europe")

    # Create City
    net.Add_Qonnector("Qonnector")
    net.Add_Qlient("QlientAlice",dist_cities/2,"Qonnector")
    net.Add_Qlient("QlientBob",dist_cities/2,"Qonnector")
    # Create channels

    # Get node instances
    Qonn = net.network.get_node("Qonnector")
    Alice = net.network.get_node("QlientAlice")
    Bob = net.network.get_node("QlientBob")

    # Entanglement based QDK from city 1 to drone 1
    send = SendEPR(Alice, Bob, EPR_succ, Qonn)
    send.start()

    receive1 = ReceiveProtocol(Qonn, Qonnector_meas_succ, Qonnector_meas_flip, True, Alice)
    receive1.start()
    receive2 = ReceiveProtocol(Qonn, Qonnector_meas_succ, Qonnector_meas_flip, True, Bob)
    receive2.start()

    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display results 
    L1 = Sifting(Alice.keylist, Bob.keylist)
    sent = len(Qonn.QlientKeys[Alice.name])
    rec = len(L1)
    eff = rec/sent
    rate = ratesources*sourceeff*eff
    skr = rate*(1-2*h(QBER))

    print("Channel length: "+str(dist_cities)+"km")
    print("EPR pairs sent from center: " + str(sent))
    print("Number of qubits received in both city: " +str(rec) )
    print("Channel efficiency : "+str(eff) + " bits per channel use")
    print("Secret key rate: "+ str(skr) + " bits per second")

    return eff,rate,skr 

def Study(dist_cities):
    res = []
    eff1,rate1,skr1 = GetRatefree(dist_cities)
    eff2,rate2,skr2 =GetRatefiber(dist_cities)
    res.append([dist_cities,rate1,rate2])
    return res

mlp.set_start_method('fork')
pool_threads = os.cpu_count() - 1
pool = mlp.Pool(pool_threads)
trans = pool.map(fnct.partial(Study), dist_cities)
pool.close()
pool.join() 

## Data saving 

save_path = '../data/'    

NameFree = os.path.join(save_path, "EPRFreespace.txt")
NameFiber = os.path.join(save_path, "EPRFiber.txt")
Free01 = open(NameFree,"w")
Fiber01 = open(NameFiber,"w")


for height in trans:
    for rx in height:
        Free01.write(str(rx[1])+ "\n") 
        Fiber01.write(str(rx[2])+ "\n")

Free01.close()
Fiber01.close()