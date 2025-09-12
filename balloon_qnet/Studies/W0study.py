from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import DownlinkChannel, compute_channel_length, CachedChannel
import multiprocessing as mlp
import os
import functools as fnct
from pathlib import Path

'''This script calculates the mean transmittance of a Balloon-to-Ground vertical downlink channel for different altitudes of the balloon. 
    It creates 4 W0Theo0X.txt files with the theoretical mean transmittance of the channel and 4 W0Simu0X.txt files with 
    the simulated mean transmittance, for each of the considered initial beam waists.
'''


# Parameters
wavelength = 1550e-9
zenith_angle = 0
ground_station_alt = 0.020 #Altitude of the receiving telescope
W0 = 0.2 #Initial Beam Waist 
obs_ratio_ground = 0.3 #Obscuration ratio of the receiving telescope 
n_max_ground = 6 #Maximum radial index of correction of AO system 
Cn0 = 9.6*10**(-14) #Reference index of refraction structure constant at ground level
u_rms = 10 #Wind speed
pointing_error = 1e-6 #Pointing error variance
Qonnector_meas_succ = 0.85 #Detector efficiency at the receiver
tracking_efficiency = 0.8 #Tracking efficiency
h_balloons = range(18,38) #Altitude range of the balloon
rx_aperture_ground = 0.4 #Aperture of the receiving telescope

W0 = [0.05,0.1,0.15,0.2]#Initial Beam Waist range
simtime = 500


#Theoretical mean transmittance 

def heightTheo(h_balloons,W0):
    
    height_balloon = compute_channel_length(ground_station_alt, h_balloons, zenith_angle)
    transmittance_down = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    
    downlink_channel = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, h_balloons, 
                                       zenith_angle, pointing_error,tracking_efficiency= tracking_efficiency,Tatm=transmittance_down)
    
    eta = np.arange(1e-7, 1, 0.001)
    mean = downlink_channel._compute_mean_channel_efficiency(eta, height_balloon, detector_efficiency = Qonnector_meas_succ)
    print("Theoretical Mean: "+ str(mean)) 
    return mean
    
    
#Simulated mean transmittance

def heightSimu(h_balloons,W0):
    
    height_balloon = compute_channel_length(ground_station_alt, h_balloons, zenith_angle)
    transmittance_down = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    # Initialize network
    net = QEurope("Europe")

    # Create quantum City 1
    net.Add_Qonnector("QonnectorCity1")

    # Create drone 1
    net.Add_Qonnector("QonnectorDroneCity1")

    # Create channels
    downlink_channel = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n_max_ground, Cn0, u_rms, wavelength, ground_station_alt, h_balloons, zenith_angle, 
                                       pointing_error,tracking_efficiency= tracking_efficiency,Tatm=transmittance_down)
    a = downlink_channel._compute_loss_probability(height_balloon,math.ceil(simtime/Qonnector_init_time ))

    downlink = CachedChannel(a)
    
    # Connect nodes
    net.connect_qonnectors("QonnectorCity1", "QonnectorDroneCity1", distance = height_balloon, loss_model = downlink)

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
    print("W0 : " + str(W0) + "m")
    print("Number of qubits sent by the Balloon: " +str(len(balloon.QlientKeys[city.name])) )
    print("Number of qubits received by Bob (Qlient): " +str(len(city.QlientKeys[balloon.name])) )
    print("Channel efficiency : "+str(chan_eff) + " bits per channel use")
    
    return chan_eff

    
#Function to parallelize

def Study(h_balloons):
    res = []
    for w in W0:
        res.append([h_balloons,w,'simu',heightSimu(h_balloons,w),'theo',heightTheo(h_balloons,w)])
    return res

#Parallelized version
mlp.set_start_method('fork')
pool_threads = os.cpu_count() - 1
pool = mlp.Pool(pool_threads)
trans = pool.map(fnct.partial(Study), h_balloons)
pool.close()
pool.join() 

#Data saving

save_path = '../data/'    

NameSimu01 = os.path.join(save_path, "W0Simu01.txt") 
NameTheo01 = os.path.join(save_path, "W0Theo01.txt")
NameSimu02 = os.path.join(save_path, "W0Simu02.txt")
NameTheo02 = os.path.join(save_path, "W0Theo02.txt")
NameSimu03 = os.path.join(save_path, "W0Simu03.txt")
NameTheo03 = os.path.join(save_path, "W0Theo03.txt")
NameSimu04 = os.path.join(save_path, "W0Simu04.txt")
NameTheo04 = os.path.join(save_path, "W0Theo04.txt")

Simu01 = open(NameSimu01, "w")
Theo01 = open(NameTheo01, "w")
Simu02 = open(NameSimu02, "w")
Theo02 = open(NameTheo02, "w")
Simu03 = open(NameSimu03, "w")
Theo03 = open(NameTheo03, "w")
Simu04 = open(NameSimu04, "w")
Theo04 = open(NameTheo04, "w")

for height in trans:
    for rx in height:
        if rx[1]==0.05:
            Simu01.write(str(rx[3])+ "\n")
            Theo01.write(str(rx[5])+ "\n")
        if rx[1]==0.1:
            Simu02.write(str(rx[3])+ "\n") 
            Theo02.write(str(rx[5])+ "\n")
        if rx[1]==0.15:
            Simu03.write(str(rx[3])+ "\n") 
            Theo03.write(str(rx[5])+ "\n")
        if rx[1]==0.2:
            Simu04.write(str(rx[3])+ "\n") 
            Theo04.write(str(rx[5])+ "\n")

Simu01.close()
Simu02.close()
Simu03.close()
Simu04.close()
Theo01.close()
Theo02.close()
Theo03.close()
Theo04.close()
