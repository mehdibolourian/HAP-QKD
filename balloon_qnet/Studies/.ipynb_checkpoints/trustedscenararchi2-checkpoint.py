from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
import balloon_qnet.cn2 as cn2
from balloon_qnet.free_space_losses import DownlinkChannel, CachedChannel, compute_height_min_horiz, HorizontalChannel, RE

"""This script estimates and print the BB84 secret key rate between the two Qlients of the Italian network,
 as described in scenario of Fig17(b) of the paper """

# Parameters
wavelength = 1550e-9
ground_station_alt = 0.020 #Altitude of the receiving telescopes
W0 = 0.1 #Initial Beam Waist 
W0_balloon = 0.2 #Initial Beam Waist 
obs_ratio_ground = 0.3 #Obscuration ratio of the receiving telescope on the ground
obs_ratio_drone = 0.3 #Obscuration ratio of the receiving telescope in the balloon
n_max_ground = 6 #Maximum radial index of correction of AO system 
Cn0 = 9.6*10**(-14) #Reference index of refraction structure constant at ground level
u_rms = 10 #Wind speed
pointing_error = 1e-6 #Pointing error variance
Qonnector_meas_succ = 0.85 #Detector efficiency at the receiver
drone_meas_succ = 0.25 #Detector efficiency in the balloon
tracking_efficiency = 0.8 #Tracking efficiency
h_balloons = 35 #Altitude of the balloon (km)
rx_aperture = 0.4 #Aperture of the receiving telescope (cm)

dist_cities = 189 #distance between the cities (km)
distAlice = 36 #distance between Alice and her Qonnector (km)
distBob = 50 #distance between Bob and her Qonnector (km)
zenith_angle = 0

#Parameters and function to calculate the secret key rate
ratesources = 80e6
sourceeff=0.01
QBER = 0.04

def h(p):
    """Binary entropy function"""
    return -p*np.log2(p)-(1-p)*np.log2(1-p)

simtime = 500000

def GetRate(h_balloons,dist_cities):
    ns.sim_reset()
    alpha = dist_cities/(2*RE)
    chan_length = 2*(h_balloons + RE)*np.sin(alpha)
    transmittance_vert = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    print("channel length: "+ str(chan_length))

    hmin = compute_height_min_horiz(chan_length,h_balloons)
    print("hmin : " +str(hmin))
    transmittance_horiz = transmittance.horizontal(hmin, chan_length, wavelength*1e9)

    Cn2_drone_to_drone = cn2.hufnagel_valley(hmin*10**3, u_rms, Cn0)
    # Initialize network
    net = QEurope("Europe")

    # Create two quantum Cities
    net.Add_Qonnector("QonnectorPadova")
    net.Add_Qlient("QlientAlice",distAlice,"QonnectorPadova")

    net.Add_Qonnector("QonnectorFlorence")
    net.Add_Qlient("QlientBob",distBob,"QonnectorFlorence")
    
    # Create drones
    net.Add_Qonnector("QonnectorDrone1")
    net.Add_Qonnector("QonnectorDrone2")
    # Create channels
    downlink_channel1 = DownlinkChannel(W0, rx_aperture, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, 
                                                          h_balloons, zenith_angle, pointing_error = pointing_error, tracking_efficiency = tracking_efficiency,
                                                            Tatm = transmittance_vert)
    downlink_channel2 = DownlinkChannel(W0, rx_aperture, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, 
                                                          h_balloons, zenith_angle, pointing_error = pointing_error, tracking_efficiency = tracking_efficiency,
                                                            Tatm = transmittance_vert)
    
    a = downlink_channel1._compute_loss_probability(h_balloons-ground_station_alt,math.ceil(simtime/Qonnector_init_time ))
    b = downlink_channel2._compute_loss_probability(h_balloons-ground_station_alt,math.ceil(simtime/Qonnector_init_time ))
    downlink1 = CachedChannel(a)
    downlink2 = CachedChannel(b)

    horizontal_channel = HorizontalChannel(W0_balloon, rx_aperture_drone, obs_ratio_drone, Cn2_drone_to_drone, wavelength, pointing_error, tracking_efficiency,transmittance_horiz)
    c = horizontal_channel._compute_loss_probability(chan_length,math.ceil(simtime/Qonnector_init_time ))
    horizlink = CachedChannel(c)
    # Connect nodes
    net.connect_qonnectors("QonnectorPadova", "QonnectorDrone1", distance = h_balloons-ground_station_alt, loss_model = downlink1)
    net.connect_qonnectors("QonnectorFlorence", "QonnectorDrone2", distance = h_balloons-ground_station_alt, loss_model = downlink2)
    net.connect_qonnectors("QonnectorDrone1", "QonnectorDrone2", distance = chan_length, loss_model = horizlink)

    # Get node instances
    city1 = net.network.get_node("QonnectorPadova")
    city2 = net.network.get_node("QonnectorFlorence")
    balloon1 = net.network.get_node("QonnectorDrone1")
    balloon2 = net.network.get_node("QonnectorDrone2")
    Alice = net.network.get_node("QlientAlice")
    Bob= net.network.get_node("QlientBob")

    # BB84 between each sublink
    send1 = SendBB84(city1, 1, 0, Alice)
    send1.start()
    receive1 = ReceiveProtocol(Alice, Qonnector_meas_succ, Qonnector_meas_flip, True, city1)
    receive1.start()

    send2 = SendBB84(city1, 1, 0, balloon1)
    send2.start()
    receive2 = ReceiveProtocol(balloon1, Qonnector_meas_succ, Qonnector_meas_flip, True, city1)
    receive2.start()
    
    send3 = SendBB84(city2, 1, 0, balloon2)
    send3.start()
    receive3 = ReceiveProtocol(balloon2, Qonnector_meas_succ, Qonnector_meas_flip, True, city2)
    receive3.start()
    
    send4 = SendBB84(balloon2, 1, Qonnector_init_flip, balloon1)
    send4.start()
    receive4 = ReceiveProtocol(balloon1, drone_meas_succ, Qonnector_meas_flip, True, balloon2)
    receive4.start()

    send5 = SendBB84(city2, 1, 0, Bob)
    send5.start()
    receive5 = ReceiveProtocol(Bob, Qonnector_meas_succ, Qonnector_meas_flip, True, city2)
    receive5.start()

    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display results 
    sentAlice = len(Alice.keylist)
    recfromAlice = len(city1.QlientKeys[Alice.name])
    eff1 = recfromAlice/sentAlice
    rate1 = ratesources*sourceeff*eff1
    skr1 = rate1*(1-2*h(QBER))/1000

    sentballoon1 = len(balloon1.QlientKeys[city1.name])
    recfromballoon1 = len(city1.QlientKeys[balloon1.name])
    eff2 = recfromballoon1/sentballoon1
    rate2 = ratesources*sourceeff*eff2
    skr2 = rate2*(1-2*h(QBER))/1000

    sentballoon2 = len(balloon2.QlientKeys[city2.name])
    recfromballoon2 = len(city2.QlientKeys[balloon2.name])
    eff3 = recfromballoon2/sentballoon2
    rate3 = ratesources*sourceeff*eff3
    skr3 = rate3*(1-2*h(QBER))/1000

    sentballoon = len(balloon1.QlientKeys[balloon2.name])
    recfromballoon = len(balloon2.QlientKeys[balloon1.name])
    eff4 = recfromballoon/sentballoon
    rate4 = ratesources*sourceeff*eff4
    skr4 = rate4*(1-2*h(QBER))/1000

    sentBob = len(Bob.keylist)
    recfromBob = len(city2.QlientKeys[Bob.name])
    eff5 = recfromBob/sentBob
    rate5 = ratesources*sourceeff*eff5
    skr5 = rate5*(1-2*h(QBER))/1000

    print("Alice -> Qonn efficiency "+str(eff1))
    print("Alice -> Qonn SKR "+str(skr1) +" kbit/sec")
    print("balloon1 -> Qonn1 efficiency " + str(eff2) )
    print("Balloon1 -> Qonn1 SKR "+str(skr2) +" kbit/sec")
    print("balloon2 -> Qonn2 efficiency " + str(eff3) )
    print("Balloon2 -> Qonn2 SKR "+str(skr3) +" kbit/sec")
    print("balloon1 -> Balloon2 efficiency " + str(eff4) )
    print("Balloon1 -> Balloon2 SKR "+str(skr4) +" kbit/sec")
    print("Bob -> Qonn efficiency "+str(eff5))
    print("Bob -> Qonn SKR "+str(skr5) +" kbit/sec")

    skrs = {"AliceQonnSKR": skr1,"balloonQonn1SKR":skr2,"balloonQonn2SKR":skr3,"balloonballoonSKR":skr4,"BobQonnSKR":skr5}
    effs={"AliceQonneff": eff1,"balloonQonn1eff":eff2,"balloonQonn2eff":eff3,"balloonballooneff":eff4,"BobQonneff":eff5}
    return skrs 

def AverageRates(n):
    tot1 = []
    tot2 = []
    tot3 = []
    tot4 = []
    tot5 = []
    for i in range(n):
        skr = GetRate(h_balloons,dist_cities)
        tot1.append(skr["AliceQonnSKR"])
        tot2.append(skr["balloonQonn1SKR"])
        tot3.append(skr["balloonQonn2SKR"])
        tot4.append(skr["balloonballoonSKR"])
        tot5.append(skr["BobQonnSKR"])
    return tot1, tot2, tot3, tot4, tot5

tot1, tot2, tot3, tot4, tot5 = AverageRates(20)  
rateAlicetoQonn = np.mean(tot1)
rateAlicetoQonnstd = np.std(tot1)
rateBalloontoQonn1 = np.mean(tot2)
rateBalloontoQonn1std = np.std(tot2)
rateBalloontoQonn2 = np.mean(tot3)
rateBalloontoQonn2std = np.std(tot3)
rateBalloontoballoon = np.mean(tot4)
rateBalloontoballoonstd = np.std(tot4)
rateBobtoQonn = np.mean(tot5)
rateBobtoQonnstd = np.std(tot5)
print("BB84 rates  \n")
print("Alice to Qonn (fiber): \n")
print("Rate: " + str(rateAlicetoQonn) + "+/-" + str(rateAlicetoQonnstd) + "\n")
print("Bob to Qonn (fiber): \n")
print("Rate: " + str(rateBobtoQonn) + "+/-" + str(rateBobtoQonnstd) + "\n")
print("Balloon to Qonn1 (free space)\n")
print("Rate: " + str(rateBalloontoQonn1) + "+/-" + str(rateBalloontoQonn1std) + "\n")
print("Balloon to Qonn2 (free space)\n")
print("Rate: " + str(rateBalloontoQonn2) + "+/-" + str(rateBalloontoQonn2std) + "\n")
print("Balloon to balloon (horizontal free space)\n")
print("Rate: " + str(rateBalloontoballoon) + "+/-" + str(rateBalloontoballoonstd) + "\n")