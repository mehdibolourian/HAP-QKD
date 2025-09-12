from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import UplinkChannel, CachedChannel, RE
import multiprocessing as mlp
import os
import functools as fnct

# Parameters
wavelength = 1550e-9
EPR_succ = 1
ground_station_alt = 0.020 #Altitude of the receiving telescopes
W0 = 0.1 #Initial Beam Waist 
obs_ratio_ground = 0.3 #Obscuration ratio of the receiving telescope same
n_max_ground = 10 #Maximum radial index of correction of AO system 
Cn0 = 9.6*10**(-14) #Reference index of refraction structure constant at ground level
u_rms = 10 #Wind speed
pointing_error = 1e-6 #Pointing error variance
Qonnector_meas_succ = 0.85 #Detector efficiency at the receiver
tracking_efficiency = 0.8 #
h_balloons = 35 #Altitude of the balloon (km)
rx_aperture = 0.4 #Aperture of the receiving telescope (cm)
p_transmit=1

dist_cities = [10,20,]#30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200] #distance between the cities (km)
ratesources = 80e6
sourceeff=0.01
QBER = 0.04

def h(p):
    """Binary entropy function"""
    return -p*np.log2(p)-(1-p)*np.log2(1-p)


simtime = 100

def GetRateMDI(dist_cities):
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
    uplink_channel1 = UplinkChannel(W0, rx_aperture, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, 
                                                          h_balloons, zenith_angle, pointing_error = pointing_error, tracking_efficiency = tracking_efficiency,
                                                            Tatm = transmittance_vert)
    uplink_channel2 = UplinkChannel(W0, rx_aperture, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, 
                                                          h_balloons, zenith_angle, pointing_error = pointing_error, tracking_efficiency = tracking_efficiency,
                                                            Tatm = transmittance_vert)
    
    a = uplink_channel1._compute_loss_probability(chan_length,math.ceil(simtime/Qonnector_init_time ))
    b = uplink_channel2._compute_loss_probability(chan_length,math.ceil(simtime/Qonnector_init_time ))
    uplink1 = CachedChannel(a)
    uplink2 = CachedChannel(b)
    # Connect nodes
    net.connect_qonnectors( "QonnectorDrone","QonnectorPadova", distance = chan_length, loss_model = uplink1)
    net.connect_qonnectors("QonnectorFlorence", "QonnectorDrone", distance = chan_length, loss_model = uplink2)

    # Get node instances
    city1 = net.network.get_node("QonnectorPadova")
    city2 = net.network.get_node("QonnectorFlorence")
    balloon = net.network.get_node("QonnectorDrone")

    # Entanglement based QDK from city 1 to drone 1
    send1 = SendBB84(balloon,1,0,city1)
    send1.start()
    send2 = SendBB84(balloon,1,0,city2)
    send2.start()

    receive1 = BSMProtocol(city1, city2,0.5*0.25*0.25,balloon)
    receive1.start()

    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display results 
    rec = len(balloon.QlientKeys[city1.name])
    sent = len(city1.QlientKeys[balloon.name])
    eff = rec/sent
    rate = ratesources*sourceeff*eff
    skr = rate*(1-2*h(QBER))

    print("Channel length: "+str(chan_length)+"km")
    print("Zenith angle: " + str(zenith_angle) + ' degree')
    print("BB84 state sent from city1: " + str(sent))
    print("BB84 state sent from city2: " + str(len(city2.QlientKeys[balloon.name])))
    print("Number of successfull MDI round: " +str(rec) )
    print("Channel efficiency : "+str(eff) + " bits per channel use")
    print("Secret key rate: "+ str(skr) + " bits per second")

    return eff,rate, skr


def GetRateMDIfiber(dist_cities):
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
    send1 = SendBB84(Qonn,1,0,Alice)
    send1.start()
    send2 = SendBB84(Qonn,1,0,Bob)
    send2.start()

    receive1 = BSMProtocol(Alice, Bob,0.5*0.25*0.25,Qonn)
    receive1.start()

    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display results 
    rec = len(Qonn.QlientKeys[Alice.name])
    sent = len(Alice.keylist)
    eff = rec/sent
    rate = ratesources*sourceeff*eff
    skr = rate*(1-2*h(QBER))

    print("Channel length: "+str(dist_cities)+"km")
    print("BB84 sent from Alice: " + str(sent))
    print("Number of successfull MDI round: " +str(rec) )
    print("Channel efficiency : "+str(eff) + " bits per channel use")
    print("Secret key rate: "+ str(skr) + " bits per second")

    return eff,rate,skr 

def Study(dist_cities):
    res = []
    eff1,rate1,skr1 = GetRateMDI(dist_cities)
    eff2,rate2,skr2 = GetRateMDIfiber(dist_cities)
    res.append([dist_cities,rate1,rate2])
    return res

mlp.set_start_method('fork')
pool_threads = os.cpu_count() - 1
pool = mlp.Pool(pool_threads)
trans = pool.map(fnct.partial(Study), dist_cities)
pool.close()
pool.join() 

#Data saving
save_path = '../data/'    

NameFree = os.path.join(save_path, "MDIfree.txt")
NameFiber = os.path.join(save_path, "MDIfiber.txt")
Free01 = open(NameFree,"w")
Fiber01 = open(NameFiber,"w")

for height in trans:
    for rx in height:
        Fiber01.write(str(rx[2])+ "\n") 
        Free01.write(str(rx[1])+ "\n")

Free01.close()
Fiber01.close()
