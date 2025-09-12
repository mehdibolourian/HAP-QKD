from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
from balloon_qnet.free_space_losses import DownlinkChannel, compute_channel_length, CachedChannel
import multiprocessing as mlp
import os
import functools as fnct

'''This script calculates the mean transmittance of a Balloon-to-Ground vertical downlink channel for different altitudes of the balloon and 
   different order of correction of the AO system. 
   It creates 9 AOballoonTheoX.txt files with the theoretical mean transmittance of the channel and 9 AOballoonTheoX.txt files with 
   the simulated mean transmittance, where X ∈ [2,3,4,5,6,7,8,9,10] is the maximum radial index of correction of AO system .
'''

# Parameters
wavelength = 1550e-9
zenith_angle = 0
ground_station_alt = 0.020 #Altitude of the receiving telescope
obs_ratio_ground = 0.3 #Obscuration ratio of the receiving telescope same
Cn0 = 9.6*10**(-14) #Reference index of refraction structure constant at ground level
u_rms = 10 #Wind speed
pointing_error = 1e-6 #Pointing error variance
Qonnector_meas_succ = 0.85 #Detector efficiency at the receiver
tracking_efficiency = 0.8 # Tracking efficiency
rx_aperture_ground = 0.3 #Aperture of the receiving telescope 
W0 = 0.1 #Initial Beam Waist 

h_balloons = range(18,38) #Altitude range of the balloon
n_max_ground = [2,3,4,5,6,7,8,9,10] #Maximum radial index of correction of AO system 

simtime = 500000 #Simulation time

#Theoretical mean transmittance 

def heightTheo(h_balloons,n):
    
    height_balloon = compute_channel_length(ground_station_alt, h_balloons, zenith_angle)
    transmittance_down = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    
  
    downlink_channel = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n, Cn0, u_rms, wavelength, ground_station_alt, h_balloons, zenith_angle,
                                            pointing_error = pointing_error, tracking_efficiency = tracking_efficiency, Tatm = transmittance_down)
    eta = np.arange(1e-7, 1, 0.001)
    mean = downlink_channel._compute_mean_channel_efficiency(eta, height_balloon, detector_efficiency = Qonnector_meas_succ)
    print("Theoretical mean " + str(mean) + " with AO " + str(n))
    return mean
    
    
#Simulated mean transmittance

def heightSimu(h_balloons,n):
    
    height_balloon = compute_channel_length(ground_station_alt, h_balloons, zenith_angle)
    transmittance_down = transmittance.slant(ground_station_alt, h_balloons, wavelength*1e9, zenith_angle)
    # Initialize network
    net = QEurope("Europe")

    # Create quantum City 1
    net.Add_Qonnector("QonnectorCity1")

    # Create drone 1
    net.Add_Qonnector("QonnectorDroneCity1")

    downlink_channel = DownlinkChannel(W0, rx_aperture_ground, obs_ratio_ground, n, Cn0, u_rms, wavelength, ground_station_alt, h_balloons, zenith_angle,
                                            pointing_error = pointing_error, tracking_efficiency = tracking_efficiency, Tatm = transmittance_down)
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

    print("Height of the balloon: "+str(h_balloons)+"km")
    print("AO : " + str(n))
    print("Number of qubits sent by the Balloon: " +str(len(balloon.QlientKeys[city.name])) )
    print("Number of qubits received by Bob (Qlient): " +str(len(city.QlientKeys[balloon.name])) )
    print("Channel efficiency : "+str(chan_eff) + " bits per channel use")
    
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

#Data saving    
save_path = '../data/'    

NameSimu02 = os.path.join(save_path, "AOballoonSimu2.txt")
NameTheo02 = os.path.join(save_path, "AOballoonTheo2.txt")

NameSimu03 = os.path.join(save_path, "AOballoonSimu3.txt")
NameTheo03 = os.path.join(save_path, "AOballoonTheo3.txt")

NameSimu04 = os.path.join(save_path, "AOballoonSimu4.txt")
NameTheo04 = os.path.join(save_path, "AOballoonTheo4.txt")

NameSimu05 = os.path.join(save_path, "AOballoonSimu5.txt")
NameTheo05 = os.path.join(save_path, "AOballoonTheo5.txt")

NameSimu06 = os.path.join(save_path, "AOballoonSimu6.txt")
NameTheo06 = os.path.join(save_path, "AOballoonTheo6.txt")

NameSimu07 = os.path.join(save_path, "AOballoonSimu7.txt")
NameTheo07 = os.path.join(save_path, "AOballoonTheo7.txt")

NameSimu08 = os.path.join(save_path, "AOballoonSimu8.txt")
NameTheo08 = os.path.join(save_path, "AOballoonTheo8.txt")

NameSimu09 = os.path.join(save_path, "AOballoonSimu9.txt")
NameTheo09 = os.path.join(save_path, "AOballoonTheo9.txt")

NameSimu10 = os.path.join(save_path, "AOballoonSimu10.txt")
NameTheo10 = os.path.join(save_path, "AOballoonTheo10.txt")


Simu02 = open(NameSimu02, "w")
Theo02 = open(NameTheo02, "w")
Simu03 = open(NameSimu03, "w")
Theo03 = open(NameTheo03, "w")
Simu04 = open(NameSimu04, "w")
Theo04 = open(NameTheo04, "w")
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
        if rx[1]==2:
            Simu02.write(str(rx[3])+ "\n") 
            Theo02.write(str(rx[5])+ "\n")
        if rx[1]==3:
           Simu03.write(str(rx[3])+ "\n") 
           Theo03.write(str(rx[5])+ "\n")
        if rx[1]==4:
           Simu04.write(str(rx[3])+ "\n") 
           Theo04.write(str(rx[5])+ "\n")
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


Simu02.close()
Simu03.close()
Simu04.close()
Simu05.close()
Simu06.close()
Simu07.close()
Simu08.close()
Simu09.close()
Simu10.close()

Theo02.close()
Theo03.close()
Theo04.close()
Theo05.close()
Theo06.close()
Theo07.close()
Theo08.close()
Theo09.close()
Theo10.close()