from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import UplinkChannel, compute_channel_length, CachedChannel
import multiprocessing as mlp
import os
import functools as fnct

'''This script calculates the mean transmittance of a Ground-To-Balloon vertical uplink channel for different altitudes of the balloon and 
   different order of correction of the AO system. 
   It creates 6 UplinkAOTheoX.txt files with the theoretical mean transmittance of the channel and 6 UplinkAOSimuX.txt files with 
   the simulated mean transmittance, where X ∈ [5,6,7,8,9,10] is the maximum radial index of correction of AO system .
'''

# Parameters
wavelength = 1550e-9
zenith_angle = 0
ground_station_alt = 0.020 #Altitude of the receiving telescope
W0 = 0.1 #Initial Beam Waist
obs_ratio_ground = 0.3 #Obscuration ratio of the receiving telescope 
Cn0 = 9.6*10**(-14) #Reference index of refraction structure constant at ground level
u_rms = 10 #Wind speed
Qonnector_meas_succ = 0.25 #Detector efficiency at the receiver
tracking_efficiency=0.8 #Tracking efficiency
h_balloons = range(18,38) #Altitude range of the balloon
tx_aperture= 0.3 #Aperture of the receiving telescope in the balloon
pointing_error = 1e-6 #Pointing error variance

n_max_ground = [5,6,7,8,9,10] #range of maximal radial index correction of the AO system
simtime = 500000 #Simulation time

#Theoretical mean transmittance 

def heightTheo(h_balloons, n_max_ground):
    
    height_balloon = compute_channel_length(ground_station_alt, h_balloons, zenith_angle)
    transmittance_up = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    
    uplink_channel = UplinkChannel(W0, tx_aperture, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, h_balloons, zenith_angle, 
                                   pointing_error = pointing_error, tracking_efficiency = tracking_efficiency, Tatm = transmittance_up)
    
    eta = np.arange(1e-7, 1, 0.001)
    mean = uplink_channel._compute_mean_channel_efficiency(eta, height_balloon, detector_efficiency = Qonnector_meas_succ)
    print("Theoretical Mean: "+ str(mean)) 
    return mean
    
    
#Simulated mean transmittance

def heightSimu(h_balloons,n_max_ground):
    
    height_balloon = compute_channel_length(ground_station_alt, h_balloons, zenith_angle)
    transmittance_up = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    # Initialize network
    net = QEurope("Europe")

    # Create quantum City 1
    net.Add_Qonnector("QonnectorCity1")

    # Create drone 1
    net.Add_Qonnector("QonnectorDroneCity1")

    # Create channels
    uplink_channel = UplinkChannel(W0, tx_aperture, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, h_balloons, zenith_angle, 
                                   pointing_error = pointing_error, tracking_efficiency = tracking_efficiency, Tatm = transmittance_up)
    
    a = uplink_channel._compute_loss_probability(height_balloon,math.ceil(simtime/Qonnector_init_time ))

    uplink = CachedChannel(a)
    
    # Connect nodes
    net.connect_qonnectors("QonnectorCity1", "QonnectorDroneCity1", distance = height_balloon, loss_model = uplink)

    # Get node instances
    city = net.network.get_node("QonnectorCity1")
    balloon = net.network.get_node("QonnectorDroneCity1")

    # BB84 from city 1 to drone 1
    send = SendBB84(city, Qonnector_init_succ, Qonnector_init_flip, balloon)
    send.start()
    receive = ReceiveProtocol(balloon, Qonnector_meas_succ, Qonnector_meas_flip, True, city)
    receive.start()


    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display and results 
    L1 = Sifting(balloon.QlientKeys[city.name], city.QlientKeys[balloon.name])
    chan_eff = len(city.QlientKeys[balloon.name])/len(balloon.QlientKeys[city.name])

    print("Height of the balloon : " + str(h_balloons) + "km")
    print("Aperture of the receiving telescope : " + str(tx_aperture) + "cm")
    print("Number of qubits sent by the Balloon: " +str(len(balloon.QlientKeys[city.name])) )
    print("Number of qubits received by Bob (Qlient): " +str(len(city.QlientKeys[balloon.name])) )
    print("Channel efficiency : "+str(chan_eff) + " bits per channel use")
    print("QBER : "+str(estimQBER(L1)))
    
    return chan_eff

    
#Function to parallelize

def Study(h_balloons):
    res = []
    for n in n_max_ground:
        res.append([h_balloons,n,'simu',heightSimu(h_balloons,n),'theo',heightTheo(h_balloons,n)])
    return res


#Parallelized version
mlp.set_start_method('fork')
pool_threads = os.cpu_count() - 1
pool = mlp.Pool(pool_threads)
trans = pool.map(fnct.partial(Study), h_balloons)
pool.close()
pool.join() 

#Data Saving
save_path = '../data/'    

NameSimu05 = os.path.join(save_path, "UplinkAOSimu5.txt")
NameTheo05 = os.path.join(save_path, "UplinkAOTheo5.txt")

NameSimu06 = os.path.join(save_path, "UplinkAOSimu6.txt")
NameTheo06 = os.path.join(save_path, "UplinkAOTheo6.txt")

NameSimu07 = os.path.join(save_path, "UplinkAOSimu7.txt")
NameTheo07 = os.path.join(save_path, "UplinkAOTheo7.txt")

NameSimu08 = os.path.join(save_path, "UplinkAOSimu8.txt")
NameTheo08 = os.path.join(save_path, "UplinkAOTheo8.txt")

NameSimu09 = os.path.join(save_path, "UplinkAOSimu9.txt")
NameTheo09 = os.path.join(save_path, "UplinkAOTheo9.txt")

NameSimu10 = os.path.join(save_path, "UplinkAOSimu10.txt")
NameTheo10 = os.path.join(save_path, "UplinkAOTheo10.txt")

Simu05 = open(NameSimu05, "w")
Theo05 = open(NameTheo05, "w")
Simu06 = open(NameSimu06, "w")
Theo06 = open(NameTheo06, "w")
Simu07 = open(NameSimu07, "w")
Theo07 = open(NameTheo07, "w")
Simu08 = open(NameSimu08, "w")
Theo08 = open(NameTheo08, "w")
Simu09 = open(NameSimu09, "w")
Theo09 = open(NameTheo09, "w")
Simu10 = open(NameSimu10, "w")
Theo10 = open(NameTheo10, "w")


for height in trans:
    for rx in height:
        if rx[1]==5:
           Simu05.write(str(rx[3])+ "\n") 
           Theo05.write(str(rx[5])+ "\n")
        if rx[1]==6:
           Simu06.write(str(rx[3])+ "\n") 
           Theo06.write(str(rx[5])+ "\n")
        if rx[1]==7:
           Simu07.write(str(rx[3])+ "\n") 
           Theo07.write(str(rx[5])+ "\n")
        if rx[1]==8:
           Simu08.write(str(rx[3])+ "\n") 
           Theo08.write(str(rx[5])+ "\n")
        if rx[1]==9:
           Simu09.write(str(rx[3])+ "\n") 
           Theo09.write(str(rx[5])+ "\n")
        if rx[1]==10:
           Simu10.write(str(rx[3])+ "\n") 
           Theo10.write(str(rx[5])+ "\n")


Simu05.close()
Simu06.close()
Simu07.close()
Simu08.close()
Simu09.close()
Simu10.close()

Theo05.close()
Theo06.close()
Theo07.close()
Theo08.close()
Theo09.close()
Theo10.close()