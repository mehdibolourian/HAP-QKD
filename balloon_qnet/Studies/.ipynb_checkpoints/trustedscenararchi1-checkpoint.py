from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import DownlinkChannel, compute_channel_length,CachedChannel,lut_zernike_index_pd, RE

"""This script estimates and print the BB84 secret key rate between the two Qlients of the Italian network,
 as described in scenario of Fig17(a) of the paper """

# Parameters
wavelength = 1550e-9
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

    # BB84 between each sublink
    send1 = SendBB84(city1, 1, 0, Alice)
    send1.start()
    receive1 = ReceiveProtocol(Alice, Qonnector_meas_succ, Qonnector_meas_flip, True, city1)
    receive1.start()

    send2 = SendBB84(city1, 1, 0, balloon)
    send2.start()
    receive2 = ReceiveProtocol(balloon, Qonnector_meas_succ, Qonnector_meas_flip, True, city1)
    receive2.start()
    
    send3 = SendBB84(city2, 1, 0, balloon)
    send3.start()
    receive3 = ReceiveProtocol(balloon, Qonnector_meas_succ, Qonnector_meas_flip, True, city2)
    receive3.start()
    
    send4 = SendBB84(city2, 1, 0, Bob)
    send4.start()
    receive4 = ReceiveProtocol(Bob, Qonnector_meas_succ, Qonnector_meas_flip, True, city2)
    receive4.start()

    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display results 
    sentAlice = len(Alice.keylist)
    recfromAlice = len(city1.QlientKeys[Alice.name])
    eff1 = recfromAlice/sentAlice
    rate1 = ratesources*sourceeff*eff1
    skr1 = rate1*(1-2*h(QBER))/1000

    sentballoon1 = len(balloon.QlientKeys[city1.name])
    recfromballoon1 = len(city1.QlientKeys[balloon.name])
    eff2 = recfromballoon1/sentballoon1
    rate2 = ratesources*sourceeff*eff2
    skr2 = rate2*(1-2*h(QBER))/1000

    sentballoon2 = len(balloon.QlientKeys[city2.name])
    recfromballoon2 = len(city2.QlientKeys[balloon.name])
    eff3 = recfromballoon2/sentballoon2
    rate3 = ratesources*sourceeff*eff3
    skr3 = rate3*(1-2*h(QBER))/1000

    sentBob = len(Bob.keylist)
    recfromBob = len(city2.QlientKeys[Bob.name])
    eff4 = recfromBob/sentBob
    rate4 = ratesources*sourceeff*eff4
    skr4 = rate4*(1-2*h(QBER))/1000

    print("Alice -> Qonn efficiency "+str(eff1))
    print("Alice -> Qonn SKR "+str(skr1) +" kbit/sec")
    print("balloon -> Qonn1 efficiency " + str(eff2) )
    print("Balloon -> Qonn1 SKR "+str(skr2) +" kbit/sec")
    print("balloon -> Qonn2 efficiency " + str(eff3) )
    print("Balloon -> Qonn2 SKR "+str(skr3) +" kbit/sec")
    print("Bob -> Qonn efficiency "+str(eff4))
    print("Bob -> Qonn SKR "+str(skr4) +" kbit/sec")

    skrs = {"AliceQonnSKR": skr1,"balloonQonn1SKR":skr2,"balloonQonn2SKR":skr3,"BobQonnSKR":skr4}
    effs={"AliceQonneff": eff1,"balloonQonn1eff":eff2,"balloonQonn2eff":eff3,"BobQonneff":eff4}
    return skrs 

def AverageRates(n):
    tot1 = []
    tot2 = []
    tot3 = []
    tot4 = []
    for i in range(n):
        skr = GetRate(h_balloons,dist_cities)
        tot1.append(skr["AliceQonnSKR"])
        tot2.append(skr["balloonQonn1SKR"])
        tot3.append(skr["balloonQonn2SKR"])
        tot4.append(skr["BobQonnSKR"])
    return tot1, tot2, tot3, tot4

tot1, tot2, tot3, tot4 = AverageRates(20)  
rateAlicetoQonn = np.mean(tot1)
rateAlicetoQonnstd = np.std(tot1)
rateBalloontoQonn1 = np.mean(tot2)
rateBalloontoQonn1std = np.std(tot2)
rateBalloontoQonn2 = np.mean(tot3)
rateBalloontoQonn2std = np.std(tot3)
rateBobtoQonn = np.mean(tot4)
rateBobtoQonnstd = np.std(tot4)
print("BB84 rates  \n")
print("Alice to Qonn (fiber): \n")
print("Rate: " + str(rateAlicetoQonn) + "+/-" + str(rateAlicetoQonnstd) + "\n")
print("Bob to Qonn (fiber): \n")
print("Rate: " + str(rateBobtoQonn) + "+/-" + str(rateBobtoQonnstd) + "\n")
print("Balloon to Qonn1 (free space)\n")
print("Rate: " + str(rateBalloontoQonn1) + "+/-" + str(rateBalloontoQonn1std) + "\n")
print("Balloon to Qonn2 (free space)\n")
print("Rate: " + str(rateBalloontoQonn2) + "+/-" + str(rateBalloontoQonn2std) + "\n")