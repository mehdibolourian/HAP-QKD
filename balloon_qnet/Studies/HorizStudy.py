from balloon_qnet.QEuropeFunctions import *
import balloon_qnet.transmittance as transmittance
import balloon_qnet.cn2 as cn2
from balloon_qnet.free_space_losses import HorizontalChannel, CachedChannel, compute_height_min_horiz
import multiprocessing as mlp
import os
import functools as fnct

'''This script calculates the mean transmittance of a Balloon-to-Balloon horizontal channel for different altitudes of the balloons. 
    It creates 4 HeightballoonTheo0X.txt files with the theoretical mean transmittance of the channel and 4 HeightballoonSimu0X.txt files with 
    the simulated mean transmittance, where X ∈ [18,23,28,33] is the height of the balloon.
'''

# Parameters
wavelength = 1550e-9
W0 = 0.1 #Initial Beam Waist 
obs_ratio_drone = 0.3 #Obscuration ratio of the receiving telescope
Cn0 = 9.6*10**(-14) #Reference index of refraction structure constant at ground level
u_rms = 10 #Wind speed
pointing_error = 1e-7 #Pointing error variance
Qonnector_meas_succ = 0.25 #Detector efficiency at the receiver
rx_aperture_drone = 0.3 #Aperture of the receiving telescope[0.1 -> 0.3]
tracking_efficiency = 0.80
h_balloons = [18,23,28,33] #Altitude range of the balloon

dist_balloons = [1,10,15,20,25,50,60,75,100,125,150,175,200,225,250,275,300,350,400] #Distance between balloons

simtime = 500000
#Theoretical mean transmittance 

def heightTheo(h_balloons,dist_balloons):
    
    hmin = compute_height_min_horiz(dist_balloons,h_balloons)
    print("hmin: " +str(hmin))
    transmittance_horiz = transmittance.horizontal(hmin, dist_balloons, wavelength*1e9)

    Cn2_drone_to_drone = cn2.hufnagel_valley(hmin*10**3, u_rms, Cn0)

    horizontal_channel = HorizontalChannel(W0, rx_aperture_drone, obs_ratio_drone, Cn2_drone_to_drone, wavelength, 
                                           pointing_error, tracking_efficiency,transmittance_horiz)

    eta = np.arange(1e-7, 1, 0.001)
    mean = horizontal_channel._compute_mean_channel_efficiency(eta, dist_balloons, detector_efficiency = Qonnector_meas_succ)
    print("Theoretical Mean: "+ str(mean)) 
    return mean
    
    
#Simulated mean transmittance

def heightSimu(h_balloons,dist_balloons):
    
    hmin = compute_height_min_horiz(dist_balloons,h_balloons)
    transmittance_horiz = transmittance.horizontal(hmin, dist_balloons, wavelength*1e9)

    Cn2_drone_to_drone = cn2.hufnagel_valley(hmin*10**3, u_rms, Cn0)


    # Initialize network
    net = QEurope("Europe")

    # Create quantum City 1
    net.Add_Qonnector("QonnectorBalloon1")

    # Create drone 1
    net.Add_Qonnector("QonnectorBalloon2")

    # Create channels
    horizontal_channel = HorizontalChannel(W0, rx_aperture_drone, obs_ratio_drone, Cn2_drone_to_drone, wavelength, pointing_error, tracking_efficiency,transmittance_horiz)
    a = horizontal_channel._compute_loss_probability(dist_balloons,math.ceil(simtime/Qonnector_init_time ))
    horizlink = CachedChannel(a)
    
    # Connect nodes
    net.connect_qonnectors("QonnectorBalloon1", "QonnectorBalloon2", distance = dist_balloons, loss_model = horizlink)

    # Get node instances
    balloonsend = net.network.get_node("QonnectorBalloon1")
    balloonrec = net.network.get_node("QonnectorBalloon2")

    # BB84 from city 1 to drone 1
    send = SendBB84(balloonrec, 1, Qonnector_init_flip, balloonsend)
    send.start()
    receive = ReceiveProtocol(balloonsend, Qonnector_meas_succ, Qonnector_meas_flip, True, balloonrec)
    receive.start()


    # Run simulation
    stat = ns.sim_run(duration = simtime)

    # Display and results 
    L1 = Sifting(balloonsend.QlientKeys[balloonrec.name], balloonrec.QlientKeys[balloonsend.name])
    sent = len(balloonsend.QlientKeys[balloonrec.name])
    rec = len(balloonrec.QlientKeys[balloonsend.name])
    chan_eff = rec/sent

    print("Height of the balloon : " + str(h_balloons) + "km")
    print("Distance between balloons : "+ str(dist_balloons)+ "km")
    print("Number of qubits sent by the Balloon 1: " +str(sent) )
    print("Number of qubits received by Balloon 2: " +str(rec) )
    print("Channel efficiency : "+str(chan_eff) + " bits per channel use")
    print("QBER : "+str(estimQBER(L1)))
    
    return chan_eff

    
#Function to parallelize

def Study(dist):
    res = []
    for h in h_balloons:
        res.append([dist,h,'simu',heightSimu(h,dist),'theo',heightTheo(h,dist)])
    return res

#Parallelized version
mlp.set_start_method('fork')
pool_threads = os.cpu_count() - 1
pool = mlp.Pool(pool_threads)
trans = pool.map(fnct.partial(Study), dist_balloons)
pool.close()
pool.join() 

#Data saving    
save_path = '../data/'    

NameSimu02 = os.path.join(save_path, "HorizSimu01.txt")
NameTheo02 = os.path.join(save_path, "HorizTheo01.txt")

NameSimu03 = os.path.join(save_path, "HorizSimu02.txt")
NameTheo03 = os.path.join(save_path, "HorizTheo02.txt")

NameSimu04 = os.path.join(save_path, "HorizSimu04.txt")
NameTheo04 = os.path.join(save_path, "HorizTheo04.txt")

NameSimu05 = os.path.join(save_path, "HorizSimu06.txt")
NameTheo05 = os.path.join(save_path, "HorizTheo06.txt")    

Simu01 = open(NameSimu02,"w")
Theo01 = open(NameTheo02,"w")
Simu02 = open(NameSimu03,"w")
Theo02 = open(NameTheo03,"w")
Simu04 = open(NameSimu04,"w")
Theo04 = open(NameTheo04,"w")
Simu06 = open(NameSimu05,"w")
Theo06 = open(NameTheo05,"w")

for height in trans:
    for rx in height:
        if rx[1]==18:
            Simu01.write(str(rx[3])+ "\n") 
            Theo01.write(str(rx[5])+ "\n")
        if rx[1]==23:
            Simu02.write(str(rx[3])+ "\n") 
            Theo02.write(str(rx[5])+ "\n")
        if rx[1]==28:
            Simu04.write(str(rx[3])+ "\n") 
            Theo04.write(str(rx[5])+ "\n")
        if rx[1]==33:
            Simu06.write(str(rx[3])+ "\n") 
            Theo06.write(str(rx[5])+ "\n")

Simu01.close()
Simu02.close()
Simu04.close()
Simu06.close()
Theo01.close()
Theo02.close()
Theo04.close()
Theo06.close()
