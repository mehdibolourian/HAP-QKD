#!/usr/bin/env python

#This file contains all the important functions for simulating a european quantum network 
# Netsquid has to be installed 

#Author : Raja Yehia (raja.yehia@gmail.com)

import netsquid as ns

import netsquid.components.instructions as instr
import netsquid.components.qprogram as qprog
import random 
from scipy.stats import bernoulli
import logging
import math
import numpy as np

from netsquid.components import Channel, QuantumChannel, QuantumMemory, ClassicalChannel
from netsquid.components.models.qerrormodels import FibreLossModel, DephaseNoiseModel
from netsquid.nodes import Node, DirectConnection
from netsquid.nodes.connections import Connection
from netsquid.protocols import NodeProtocol
from netsquid.components.models import DelayModel
from netsquid.components.models.delaymodels import FixedDelayModel, FibreDelayModel
from netsquid.components import QuantumMemory
from netsquid.qubits.state_sampler import StateSampler
from netsquid.components.qsource import QSource, SourceStatus
from netsquid.components.qprocessor import QuantumProcessor
from netsquid.nodes.network import Network
from netsquid.qubits import ketstates as ks
from netsquid.protocols.protocol import Signals
from netsquid.components.qprocessor import PhysicalInstruction
from netsquid.qubits import qubitapi as qapi
from netsquid.components.clock import Clock

import sys

#Qlient parameters
f_qubit_qlient = 100e6 #Qubit creation attempt frequency
Qlient_init_time = math.ceil(1e9/f_qubit_qlient) #time to create |0> in a Qlient node in ns
Qlient_init_succ = 0.008 #Probability that a a qubit creation succeed
Qlient_init_flip = 0 # probability that a qubit is flipped at its creation
Qlient_meas_succ= 0.95 #Probability that a measurement succeeds
Qlient_meas_flip = 1e-5 #Probability that the measurement outcome is flipped by the detectors 


#Qonnector parameters
Max_Qlient = 5 #Number of simultaneous link that the Qonnector can create 
f_qubit_qonn = 100e6 #Qubit creation attempt frequency in MHz
Qonnector_init_time = math.ceil(1e9/f_qubit_qonn) #time to create |0> in a Qonnector node in ns
Qonnector_init_succ = 1 #Probability that a qubit creation succeeds
Qonnector_init_flip = 0
Qonnector_meas_succ=0.85 #Probability that a measurement succeeds
Qonnector_meas_flip = 0 #Probability that the measurement outcome is flipped by the detectors 
switch_succ=0.9 #probability that transmitting a qubit from a qlient to another succeeds
BSM_succ = 0.9 #probability that a Bell state measurement of 2 qubits succeeds
EPR_succ=1 #probability that an initialisation of an EPR pair succeeds
f_EPR = 80e6 #EPR pair creation attempt frequency in MHz
EPR_time = math.ceil(1e9/f_EPR) # time to create a bell pair in a Qonnector node (ns)

#Detector noise parameters
DetectGate = 1e-10 #Detection gate window
DCRate = 100 #Dark count rate
BGRateDay = 5000 #Background Noise rate by day
BGRateNight = 500 #Background Noise rate by night

pdark = DCRate*DetectGate #Probability to get a dark count within the detection window
pbgday = BGRateDay*DetectGate*Qonnector_meas_succ #Probability to get a background click count within the detection window by day
pbgnight = BGRateNight*DetectGate*Qonnector_meas_succ #Probability to get a background click count within the detection window by night

pafterpulse = 0.05 #afterpulse probability

#Network parameter
fiber_coupling = 0.9 #Fiber coupling efficiency
fiber_loss=0.18 #Loss in fiber in dB/km
fiber_dephasing_rate = 0.02 #dephasing rate in the fiber (Hz)

#Satellite to Ground channel parameters
txDiv = 500e-6
sigmaPoint = 5e-6
rx_aperture_sat = 1
Cn2_sat = 0

#Free space channel parameter
wavelength = 1550*1e-9
W0 = wavelength/(txDiv*np.pi)
rx_aperture_drone = 0.3
rx_aperture_ground = 0.3
Cn2_drone_to_ground = 10e-16
Cn2_drone_to_drone = 10e-18
c = 299792.458 #speed of light in km/s
Tatm = 1

#Quantum operations accessible to the Qonnectors
qonnector_physical_instructions = [
    PhysicalInstruction(instr.INSTR_INIT, duration=Qonnector_init_time),
    PhysicalInstruction(instr.INSTR_H, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_X, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_Z, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_S, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_I, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_CNOT, duration=4, parallel=True),
    PhysicalInstruction(instr.INSTR_MEASURE, duration=1, parallel=True, topology=[0,1]),
    PhysicalInstruction(instr.INSTR_MEASURE_BELL, duration = 1, parallel=True),
    PhysicalInstruction(instr.INSTR_SWAP, duration = 1, parallel=True)
]

#Quantum operations accessible to the Qlient
qlient_physical_instructions = [
    PhysicalInstruction(instr.INSTR_INIT, duration=Qlient_init_time),
    PhysicalInstruction(instr.INSTR_H, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_X, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_Z, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_S, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_I, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_MEASURE, duration=1, parallel=False, topology=[0])
]


class Qlient_node(Node):
    """A Qlient node
    
    Parameters:
    name: name of the Qlient
    phys_instruction: list of physical instructions for the Qlient
    keylist: list of bits for QKD
    ports: list of two ports: one to send to Qonnector and one to receive
    """
    
    def __init__(self, name, phys_instruction, keylist=None,listports=None):
        super().__init__(name=name)
        qmem = QuantumProcessor("QlientMemory{}".format(name), num_positions=1,
                                phys_instructions=phys_instruction)
        self.qmemory = qmem
        self.keylist=keylist
        self.listports=listports
        
class Qonnector_node(Node):
    """A Qonnector node
    
    Parameters:
    QlientList: List of connected Qlients
    QlientPorts: Dictionnary of the form {Qlient: [port_to_send, port_to_receive]}
    QlientKeys : Dictionnary for QKD of the form {Qlient: [key]}
    type : Type of connector. Can be "ground", "satellite" or "drone".
    """
    
    def __init__(self, name, QlientList = None,
                  QlientPorts = None, QlientKeys = None, type = "ground"):
        super().__init__(name=name)
        self.QlientList = QlientList
        self.QlientPorts = QlientPorts
        self.QlientKeys = QlientKeys
        self.type = type
        

class QEurope():
    
    def __init__(self,name):
        """ Initialisation of a Quantum network
        
        Parameter:
        name: name of the network (str) /!\ Expected name should start with 'Qonnector' /!\
        """
        self.network = Network(name)
        self.name=name
    
    def Add_Qonnector(self, qonnectorname):
        """Method to add a Qonnector to the network
        
        Parameter :
        qonnectorname: name tof the Qonnector to add (str)
        """
        
        Qonnector = Qonnector_node(qonnectorname, QlientList = [], QlientPorts = {}, QlientKeys = {})
        self.network.add_node(Qonnector)


    def Add_Qlient(self, qlientname, distance, qonnectorto):
    
        """ Method to add a Qlient to the network. It creates a Quantum Processor at the Qonnector qonnectorto
         that is linked to the new Qlient through a fiber.
        
        Parameters :
        qlientname: name of the qlient to add (str)
        distance: distance from the Qonnector to the new node in km
        qonnectorto: Name of the Qonnector to attach the Qlient to (str)
        """
        
        network = self.network
        # Check that the Qonnector has space for the new qlient
        Qonnector = network.get_node(qonnectorto)
        if len(Qonnector.QlientList)==Max_Qlient:
            raise ValueError("You have reached the maximum Qlient capacity for this Qonnector.")
        
        #creates a qlient and adds it to the network
        Qlient = Qlient_node(qlientname,qlient_physical_instructions,keylist=[],listports=[])
        network.add_node(Qlient) 
        
        #Create quantum channels and add them to the network
        qchannel1 = QuantumChannel("QuantumChannelSto{}".format(qlientname),length=distance, delay=1,
                                   models={"quantum_loss_model": FibreLossModel(p_loss_init=1-fiber_coupling,
                                                                                p_loss_length=fiber_loss),
                                          "quantum_noise_model":DephaseNoiseModel(dephase_rate = fiber_dephasing_rate,time_independent=True)})
        qchannel2 = QuantumChannel("QuantumChannel{}toS".format(qlientname),length=distance, delay=1,
                                   models={"quantum_loss_model": FibreLossModel(p_loss_init=1-fiber_coupling,
                                                                                p_loss_length=fiber_loss),
                                           "quantum_noise_model":DephaseNoiseModel(dephase_rate = fiber_dephasing_rate,time_independent=True)})
        
        
        Qonn_send, Qlient_receive = network.add_connection(
            qonnectorto, qlientname, channel_to=qchannel1, label="quantumS{}".format(qlientname))
        Qlient_send, Qonn_receive = network.add_connection(
            qlientname, qonnectorto, channel_to=qchannel2, label="quantum{}S".format(qlientname))
    
        # Update the Qonnector's properties
        qmem = QuantumProcessor( "QonnectorMemoryTo{}".format(qlientname), num_positions=2 ,
                                phys_instructions=qonnector_physical_instructions)
        Qonnector.add_subcomponent(qmem)
        Qonnector.QlientList.append(qlientname)
        Qonnector.QlientPorts[qlientname] = [Qonn_send,Qonn_receive]
        Qonnector.QlientKeys[qlientname] = []
    
        #Update Qlient ports
        Qlient.listports = [Qlient_send, Qlient_receive]
        
        def route_qubits(msg):
            target = msg.meta.pop('internal', None)

            if isinstance(target, QuantumMemory):
                if not target.has_supercomponent(Qonnector):
                    raise ValueError("Can't internally route to a quantummemory that is not a subcomponent.")
                target.ports['qin'].tx_input(msg)
            else:
                Qonnector.ports[Qonn_send].tx_output(msg)
            
        # Connect the Qonnector's ports
        qmem.ports['qout'].bind_output_handler(route_qubits) #port to send to Qlient
        Qonnector.ports[Qonn_receive].forward_input(qmem.ports["qin"]) #port to receive from Qlient

        # Connect the Qlient's ports 
        Qlient.ports[Qlient_receive].forward_input(Qlient.qmemory.ports["qin"]) #port to receive from qonnector
        Qlient.qmemory.ports["qout"].forward_output(Qlient.ports[Qlient_send]) #port to send to qonnector
        
        #Classical channels on top of that
        cchannel1 = ClassicalChannel("ClassicalChannelSto{}".format(qlientname),length=distance, delay=1)
        cchannel2 = ClassicalChannel("ClassicalChannel{}toS".format(qlientname),length=distance, delay=1)
        
        network.add_connection(qonnectorto, qlientname, channel_to=cchannel1, 
                               label="ClassicalS{}".format(qlientname), port_name_node1="cout_{}".format(qlientname),
                               port_name_node2="cin")
        network.add_connection(qlientname, qonnectorto, channel_to=cchannel2, 
                               label="Classical{}S".format(qlientname), port_name_node1="cout",
                               port_name_node2="cin_{}".format(qlientname))
        
    def connect_qonnectors(self, qonnector_1_name, qonnector_2_name, distance, loss_model):
        
        # Create Qonnector objects
        network = self.network
        qonnector_1 = network.get_node(qonnector_1_name)
        qonnector_2 = network.get_node(qonnector_2_name)

        # Create dedicated processor at each Qonnector
        q_mem_1 = QuantumProcessor("QonnectorMemoryTo{}".format(qonnector_2_name), num_positions = 2,
                            phys_instructions = qonnector_physical_instructions)
        qonnector_1.add_subcomponent(q_mem_1)
        
        q_mem_2 = QuantumProcessor("QonnectorMemoryTo{}".format(qonnector_1_name), num_positions = 2,
                            phys_instructions=qonnector_physical_instructions)
        qonnector_2.add_subcomponent(q_mem_2)
        
        # Create biderctional quantum channel
        q_channel_1 = QuantumChannel("QonnectorChannelto{}".format(qonnector_1_name), length = distance, delay = 1,
                                    models = {"quantum_loss_model": loss_model})
        
        q_channel_2 = QuantumChannel("QonnectorChannelto{}".format(qonnector_2_name), length = distance, delay = 1,
                                    models = {"quantum_loss_model": loss_model})
        
        # Connect Qonnectors via the quantum channel
        qonnector_1_send, qonnector_2_receive = network.add_connection(
                qonnector_1, qonnector_2, channel_to = q_channel_1, label = "QonnectorChanTo{}".format(qonnector_2_name))
        qonnector_2_send, qonnector_1_receive = network.add_connection(
                qonnector_2, qonnector_1, channel_to = q_channel_2, label = "QonnectorChanTo{}".format(qonnector_1_name))

        # Update Qonnector 1 properties
        qonnector_1.QlientList.append(qonnector_2_name)
        qonnector_1.QlientPorts[qonnector_2_name] = [qonnector_1_send, qonnector_1_receive]
        qonnector_1.QlientKeys[qonnector_2_name] = []
        
        # Update Qonnector 2 properties
        qonnector_2.QlientList.append(qonnector_1_name)
        qonnector_2.QlientPorts[qonnector_1_name] = [qonnector_2_send, qonnector_2_receive]
        qonnector_2.QlientKeys[qonnector_1_name] = []

        # Update Qonnector 1 ports
        def route_qubits_1(msg):
            target = msg.meta.pop('internal', None)

            if isinstance(target, QuantumMemory):
                if not target.has_supercomponent(qonnector_1):
                    raise ValueError("Can't internally route to a quantum memory that is not a subcomponent.")
                target.ports['qin'].tx_input(msg)
            else:
                qonnector_1.ports[qonnector_1_send].tx_output(msg)
        
        q_mem_1.ports['qout'].bind_output_handler(route_qubits_1) 
        qonnector_1.ports[qonnector_1_receive].forward_input(q_mem_1.ports["qin"])   
        
        # Update Qonnector 2 ports
        def route_qubits_2(msg):
            target = msg.meta.pop('internal', None)

            if isinstance(target, QuantumMemory):
                if not target.has_supercomponent(qonnector_2):
                    raise ValueError("Can't internally route to a quantum memory that is not a subcomponent.")
                target.ports['qin'].tx_input(msg)
            else:
                qonnector_2.ports[qonnector_2_send].tx_output(msg)
    
        q_mem_2.ports['qout'].bind_output_handler(route_qubits_2) 
        qonnector_2.ports[qonnector_2_receive].forward_input(q_mem_2.ports["qin"])

        # Create biderectional classical channel
        c_channel_1 = ClassicalChannel("ClassicalChannelto{}".format(qonnector_1), length = distance, delay=1)
        c_channel_2 = ClassicalChannel("ClassicalChannelto{}".format(qonnector_2), length = distance, delay=1)
    
        # Connect Qonnectors via the classical channel
        network.add_connection(qonnector_2, qonnector_1, channel_to=c_channel_1, 
                        label = "Classicalto{}".format(qonnector_1), port_name_node1 = "cout_{}".format(qonnector_1_name),
                        port_name_node2 = "cin_{}".format(qonnector_2_name))
        network.add_connection(qonnector_1, qonnector_2, channel_to=c_channel_2, 
                        label = "Classicalto{}".format(qonnector_2), port_name_node1 = "cout_{}".format(qonnector_2_name),
                        port_name_node2 = "cin_{}".format(qonnector_1_name))
        
            
            
class ReceiveProtocol(NodeProtocol):
    
    """Protocol performed by a node to receive a state a measure it. Saves the outputs as well as the
    measurement basis and timestamps of reception in the list QlientKeys["name of the sending node"]   
    in the case of a Qonnector or in the list keylist in the case of a Qlient.
    
    Parameters:
     othernode: node from which a qubit is expected
     measurement_succ: probability that the measurement succeeds
     measurement_flip: probability that the detector flips the outcome (crosstalk)
     BB84: boolean indicating if we perform BB84 measurement (random choice of measurement basis)"""
    
    def __init__(self, othernode, measurement_succ, measurement_flip, BB84, node):
            super().__init__(node=node)
            self._othernode = othernode
            self._measurement_succ=measurement_succ
            self._BB84 = BB84
            self._measurement_flip = measurement_flip

    def run(self):
            if self.node.name[0:9] == 'Qonnector' or self.node.name[0:5]=='Drone':
                if self.node.name[0:9] == 'Qonnector':
                    mem = self.node.subcomponents["QonnectorMemoryTo{}".format(self._othernode.name)]
                    port = self.node.ports[self.node.QlientPorts[self._othernode.name][1]]
                elif self.node.name[0:5]=='Drone':
                    mem = self.node.subcomponents["DroneMemoryTo{}".format(self._othernode.name)]
                    port = self.node.ports[self.node.QlientPorts[self._othernode.name][1]]
                #print(port)
                while True:
                    yield self.await_port_input(port)
                    t = self.node.ports["cin_{}".format(self._othernode.name)].rx_input()
                    
                    b = bernoulli.rvs(self._measurement_succ)
                    
                    if b ==1 :
                        if self._BB84: #in case we perform BB84
                            base = bernoulli.rvs(0.5) #choose a random basis
                            if base < 0.5:
                                mem.execute_instruction(instr.INSTR_H, [0], physical = False)
                                base = "plusmoins"
                            else:
                                mem.execute_instruction(instr.INSTR_I, [0],physical = False)
                                base = "zeroun"
                        else:
                            base = None 
                        
                        m,_,_ = mem.execute_instruction(instr.INSTR_MEASURE,[0],output_key="M1")
                        yield self.await_program(mem,await_done=True,await_fail=True)
                        
                        flip = bernoulli.rvs(self._measurement_flip)
                        if (flip==1):
                            if m['M1'][0]==0:
                                m['M1'][0] =1
                            elif m['M1'][0]==1:
                                m['M1'][0]=0
                            
                        if m['M1'] is not None and t is not None and base is not None:
                            self.node.QlientKeys[self._othernode.name].append(([t.items[0],base],m['M1'][0]))
                            
                        elif m['M1'] is not None and t is not None:
                            self.node.QlientKeys[self._othernode.name].append((t.items,m['M1'][0]))
                            
                        elif m['M1'] is not None:
                            self.node.QlientKeys[self._othernode.name].append(m['M1'][0])
                    mem.reset()
                        
            
            else:
                mem = self.node.qmemory
                port = self.node.ports[self.node.listports[1]]
                
                while True:
                    yield self.await_port_input(port)
                    #print("qubit received")
                    #qubit, = mem.peek([0])
                    #print(mem.peek([0]))
                    t = self.node.ports["cin"].rx_input()
                    
                    b = bernoulli.rvs(self._measurement_succ)
                    #print(b)
                    if b ==1 :
                        if self._BB84: #in case we perform BB84
                            base = bernoulli.rvs(0.5) #choose a random basis
                            if base < 0.5:
                                mem.execute_instruction(instr.INSTR_H, [0], physical = False)
                                base = "plusmoins"
                            else:
                                mem.execute_instruction(instr.INSTR_I, [0],physical = False)
                                base = "zeroun"
                        else:
                            base = None 
                            
                        if not(mem.busy):
                            
                            
                            m,_,_ = mem.execute_instruction(instr.INSTR_MEASURE,[0],output_key="M1")
                            yield self.await_program(mem,await_done=True,await_fail=True)
                            #print("qubit measured")
                            
                            flip = bernoulli.rvs(self._measurement_flip)
                            if (flip==1):
                                if m['M1'][0]==0:
                                    m['M1'][0] =1
                                elif m['M1'][0]==1:
                                    m['M1'][0]=0
                                    
                            if m['M1'] is not None and t is not None and base is not None:
                                self.node.keylist.append(([t.items[0], base],m['M1'][0]))
                            
                            elif m['M1'] is not None and t is not None:
                                self.node.keylist.append((t.items,m['M1'][0]))
                        
                            elif m['M1'] is not None:
                                self.node.keylist.append(m['M1'][0])
                            
                    mem.reset()
                
                
                
class  TransmitProtocol(NodeProtocol):
    """Protocol performed by a Qonnector to transmit a qubit sent by a Qlient or a satellite to another Qlient
        
        Parameters
         Qlient_from: node from which a qubit is expected
         Qlient_to: node to which transmit the qubit received
         switch_succ: probability that the transmission succeeds"""
        
    def __init__(self, Qlient_from, Qlient_to, switch_succ, node=None, name=None):
                super().__init__(node=node, name=name)
                self._Qlient_from = Qlient_from
                self._Qlient_to = Qlient_to
                self._switch_succ=switch_succ
        
    def run(self):
            rec_mem = self.node.subcomponents["QonnectorMemoryTo{}".format(self._Qlient_from.name)]
            rec_port = self.node.ports[self.node.QlientPorts[self._Qlient_from.name][1]]
            sen_mem = self.node.subcomponents["QonnectorMemoryTo{}".format(self._Qlient_to.name)]
            #print(sen_mem)
            while True:
                rec_mem.reset()
                sen_mem.reset()
                
                yield self.await_port_input(rec_port)
                #print("qubit received at qonnector" )
                t = self.node.ports["cin_{}".format(self._Qlient_from.name)].rx_input()
                
                rec_mem.pop([0], skip_noise=True, meta_data={'internal': sen_mem})
                #print("qubit moved in qonnector's memory")               
                
                b = bernoulli.rvs(self._switch_succ)
                if b ==1 :
                    qubit, = sen_mem.pop([0])
                    self.node.ports["cout_{}".format(self._Qlient_to.name)].tx_output(t)
                    #print("qubit sent to node")
                    
                    
class SendEPR(NodeProtocol):
    """Protocol performed by a Qonnector or a satellite node to create a send EPR pair to two nodes, 
    each getting one qubit.
    
        Parameters:
         Qlient1 : Name of the first receiving node (str)
         Qlient2 : Name of the second receiving node (str)
         EPR_succ : success probability of the creation of the EPR pair
    
    """
    def __init__(self, Qlient_1, Qlient_2, EPR_succ, node = None, name = None):
        super().__init__(node=node, name=name)
        self._Qlient_1 = Qlient_1
        self._Qlient_2 = Qlient_2
        self._EPR_succ = EPR_succ
        
    def run(self):
        if self.node.name[0:9]== 'Satellite':
            mem1 = self.node.subcomponents["SatelliteMemoryTo{}".format(self._Qlient_1.name)]
            mem2 = self.node.subcomponents["SatelliteMemoryTo{}".format(self._Qlient_2.name)]
            port1 = self.node.ports[self.node.QlientPorts[self._Qlient_1.name][0]]
        else:
            mem1 = self.node.subcomponents["QonnectorMemoryTo{}".format(self._Qlient_1.name)]
            mem2 = self.node.subcomponents["QonnectorMemoryTo{}".format(self._Qlient_2.name)]
            port1 = self.node.ports[self.node.QlientPorts[self._Qlient_1.name][1]]
        state_sampler = StateSampler(qreprs=[ks.b11],
                                 probabilities=[1])

        qsource = QSource("qsource{}".format(self._Qlient_1.name+self._Qlient_2.name),
                          state_sampler=state_sampler,
                          num_ports=2,
                          timing_model=FixedDelayModel(delay=EPR_time),
                          status=SourceStatus.EXTERNAL)
        clock = Clock(name="clock",
                      start_delay=0,
                      models={"timing_model": FixedDelayModel(delay=EPR_time)})
        
        self.node.add_subcomponent(clock)
        self.node.add_subcomponent(qsource)
        clock.ports["cout"].connect(qsource.ports["trigger"])
        
        qsource.ports["qout0"].connect(mem1.ports["qin"])
        qsource.ports["qout1"].connect(mem2.ports["qin"])
        clock.start()
        
        while True:
            yield self.await_port_input(mem1.ports["qin"]) and self.await_port_input(mem2.ports["qin"])
            
            b = bernoulli.rvs(self._EPR_succ)
            if b==1:
                mem1.pop([0])
                self.node.ports["cout_{}".format(self._Qlient_1.name)].tx_output(clock.num_ticks)
                self.node.QlientKeys[self._Qlient_1.name].append(0)
                mem2.pop([0])
                self.node.ports["cout_{}".format(self._Qlient_2.name)].tx_output(clock.num_ticks)
                self.node.QlientKeys[self._Qlient_2.name].append(0)
                
            mem1.reset()
            mem2.reset()
            
class SendBB84(NodeProtocol):
    
    """Protocol performed by a node to send a random BB84 qubit |0>, |1>, |+> or |-> .
    
    Parameters:
     othernode: name of the receiving node (str).
     init_succ: probability that a qubit creation attempt succeeds.
     init_flip : probability that a qubit created is flipped before the sending.
     """
    
    def __init__(self,othernode, init_succ, init_flip,node):
        super().__init__(node=node)
        self._othernode = othernode
        self._init_succ = init_succ
        self._init_flip = init_flip
    
    def run(self):
        if self.node.name[0:9] == 'Qonnector' or self.node.name[0:9]== 'Satellite' or self.node.name[0:5]=='Drone':
            
            if self.node.name[0:9]== 'Satellite':
                mem = self.node.subcomponents["SatelliteMemoryTo{}".format(self._othernode.name)]
            elif self.node.name[0:5]=='Drone':
                mem = self.node.subcomponents["DroneMemoryTo{}".format(self._othernode.name)]
            else:
                mem = self.node.subcomponents["QonnectorMemoryTo{}".format(self._othernode.name)]
        
            
            clock = Clock(name="clock",
                      start_delay=0,
                      models={"timing_model": FixedDelayModel(delay=Qonnector_init_time )})
            self.node.add_subcomponent(clock)
            clock.start()
        
            while True:
                mem.reset()

                mem.execute_instruction(instr.INSTR_INIT,[0])
                yield self.await_program(mem,await_done=True,await_fail=True)
                #print("qubit created")
                succ = bernoulli.rvs(self._init_succ)
                if (succ == 1):                    
                    flip = bernoulli.rvs(self._init_flip)
                    if (flip == 1):
                        mem.execute_instruction(instr.INSTR_X, [0], physical = False)
            
                    base = bernoulli.rvs(0.5) #random choice of a basis
                    if base <0.5:
                        mem.execute_instruction(instr.INSTR_H,[0])
                        base = "plusmoins"
                    else:
                        mem.execute_instruction(instr.INSTR_I,[0])
                        base = "zeroun"
                
                    yield self.await_program(mem,await_done=True,await_fail=True)
                
                    t = clock.num_ticks
                    bit = bernoulli.rvs(0.5) #random choice of a bit
                    if bit < 0.5:
                        mem.execute_instruction(instr.INSTR_I, [0], physical=False)
                        self.node.QlientKeys[self._othernode.name].append(([t,base],0))
                    else:
                        if base == "zeroun":
                            mem.execute_instruction(instr.INSTR_X, [0], physical=False)
                        elif base == "plusmoins":
                            mem.execute_instruction(instr.INSTR_Z, [0], physical=False)
                        self.node.QlientKeys[self._othernode.name].append(([t,base],1))
                
                    qubit, = mem.pop([0])
                    self.node.ports["cout_{}".format(self._othernode.name)].tx_output(t)
                
                
        else:
            mem = self.node.qmemory
            clock = Clock(name="clock",start_delay=0,
                      models={"timing_model": FixedDelayModel(delay=Qlient_init_time )})
            
            self.node.add_subcomponent(clock)
            clock.start()
            
            while True:
                mem.reset()

                mem.execute_instruction(instr.INSTR_INIT,[0])
                yield self.await_program(mem,await_done=True,await_fail=True)
                    #print("qubit created")
                succ = bernoulli.rvs(self._init_succ)
                if (succ == 1):      
                    flip = bernoulli.rvs(self._init_flip)
                    if (flip == 1):
                        mem.execute_instruction(instr.INSTR_X, [0], physical = False)
            
                    base = bernoulli.rvs(0.5) #random choice of a basis
                    if base <0.5:
                        mem.execute_instruction(instr.INSTR_H,[0])
                        base = "plusmoins"
                    else:
                        mem.execute_instruction(instr.INSTR_I,[0])
                        base = "zeroun"
            
                    yield self.await_program(mem,await_done=True,await_fail=True)
                
                    t = clock.num_ticks
                    bit = bernoulli.rvs(0.5) #random choice of a bit
                    if bit < 0.5:
                        mem.execute_instruction(instr.INSTR_I, [0], physical=False)
                        self.node.keylist.append(([t,base],0))
                    else:
                        if base == "zeroun":
                            mem.execute_instruction(instr.INSTR_X, [0], physical=False)
                        elif base == "plusmoins":
                            mem.execute_instruction(instr.INSTR_Z, [0], physical=False)
                        self.node.keylist.append(([t,base],1))
            
                    qubit, = mem.pop([0])
                    self.node.ports["cout"].tx_output(t)
                
                    #print("qubit sent")

                    
                    
class BellMeasurementProgram(qprog.QuantumProgram):
    """Program to perform a Bell measurement on two qubits.

    Measurement results are stored in output keys "M1" and "M2"

    """
    default_num_qubits = 2

    def program(self):
        q1, q2 = self.get_qubit_indices(2)
        self.apply(instr.INSTR_CNOT, [q1, q2],physical=False)
        self.apply(instr.INSTR_H, [q1],physical=False)
        self.apply(instr.INSTR_MEASURE, [q1],physical=False, output_key="M1")
        self.apply(instr.INSTR_MEASURE, [q2],physical=False, output_key="M2")
        yield self.run()
        
class BSMProtocol(NodeProtocol):
    """Bell state measurement performed by a Qonnector node on two qubits received from two different nodes.
    Can only be performed by a Qonnector node. Outputs are stored in QlientKeys['name of the Qlients'].
    
    Parameters
     Qlient_1, Qlient_2: name of nodes from wich the Qonnector expects qubits (str)
     BSM_succ : Success probability of the Bell state measurement.
     """
    
    def __init__(self, Qlient_1, Qlient_2, BSM_succ, node = None, name=None):
        super().__init__(node=node, name=name)
        self._Qlient_1 = Qlient_1
        self._Qlient_2 = Qlient_2
        self._BSM_succ = BSM_succ
        
    def run(self):
        
        if not(self.node.name[0:9] == 'Qonnector'):
            return "This node cannot perform this protocol"
        
        rec_mem1 = self.node.subcomponents["QonnectorMemoryTo{}".format(self._Qlient_1.name)]
        rec_port1 = self.node.ports[self.node.QlientPorts[self._Qlient_1.name][1]]
        rec_mem2 = self.node.subcomponents["QonnectorMemoryTo{}".format(self._Qlient_2.name)]
        rec_port2 = self.node.ports[self.node.QlientPorts[self._Qlient_2.name][1]]
        
        measure_program = BellMeasurementProgram()
        
        while True:
            yield self.await_port_input(rec_port1) 
            t = self.node.ports["cin_{}".format(self._Qlient_1.name)].rx_input()
            # print("qubit1 received")
            # print(rec_mem1.peek([0]))
            # print(t)
            yield self.await_port_input(rec_port2)
            # print("qubit2 received")
            t2 = self.node.ports["cin_{}".format(self._Qlient_2.name)].rx_input()
            # print(rec_mem2.peek([0]))
            # print(t2)
            
            if t.items[0]==(t2.items[0]-1):
                
                rec_mem1.execute_instruction(instr.INSTR_SWAP, [0,1])
                yield self.await_program(rec_mem1,await_done=True,await_fail=True)
                #print("qubit 1 moved")
                #print(rec_mem1.peek([0]))
                #print(rec_mem1.peek([1]))
                
                rec_mem2.pop([0], skip_noise=True, meta_data={'internal': rec_mem1})
                yield self.await_port_input(rec_port1) 
                #print(rec_mem1.peek([0]))
                #print(rec_mem1.peek([1]))
                #print("qubit 2 moved")
                
                yield rec_mem1.execute_program(measure_program)
            
                #print(measure_program.output)
                m1, = measure_program.output["M1"]
                m2, = measure_program.output["M2"]
                #print(m1)
                #print("measurement done")
            
                b = bernoulli.rvs(self._BSM_succ)
                if m1 is not None and m2 is not None and b ==1:
                    self.node.QlientKeys[self._Qlient_1.name].append((t.items[0],m1))
                    self.node.QlientKeys[self._Qlient_2.name].append((t2.items[0]-1,m2))
            
            rec_mem1.reset()
            rec_mem2.reset() 
            
            
#Classical post-processing functions
def getTime(t):
    a,b = t
    return a[0]

def addDarkCounts(L,pdark, K):
    """Function to add dark counts to a list of outcomes. With probability pdark it will add an outcome to a 
     timestep where nothing was measured and with probability pdark/2 it will discard an outcome
     
    Parameters:
        L : List of outcomes
        pdark: probability of getting a dark count at a particular timestep
        K : Last timestep"""
    i = 0
    listestep = []
    for j in L:
        a, b = j
        listestep.append(a[0])
    while i< K:
        if i not in listestep:
            b = bernoulli.rvs(pdark)
            if b == 1:
                randbit = bernoulli.rvs(0.5)
                base = bernoulli.rvs(0.5)
                if base < 0.5:
                    L.append(([i,"plusmoins"],randbit))
                else:
                    L.append(([i,"zeroun"],randbit))
            i=i+1
        else:
            i=i+1
    L.sort(key=getTime)
    for e in L:
        if bernoulli.rvs(pdark/2)==1:
            L.remove(e)
            
            
def Sifting(Lalice, Lbob):
    """Sifting function to get a list of matching received qubit. If BB84 then the resulting list contains 
    the qubits that were sent and measured in the same basis. If EPR then the resulting list contains the qubit 
    measured by Alice and Bob that came from the same EPR pair
     
     Parameters:
     Lalice, Lbob: lists of outcomes """
    Lres = []
    for i in range(len(Lalice)):
        ta, ma = Lalice[i]
        for j in range(len(Lbob)):
            tb, mb = Lbob[j]
            if ta == tb:
                Lres.append((ma,mb))
        
    return Lres

        
def estimQBER(L):
    """Function to estimate the QBER from a list of couple (qubit sent,qubit measured)"""
    if L != []:
        Lres = []
        for i in L:
            (a,b)=i
            if a != b:
                Lres.append(b)
        return len(Lres)/len(L)
    else:
        return 1

def estimQBEREPR(L):
    """Function to estimate the QBER from a list of couple during BBM92(qubit at ALICE ,qubit at BOB)"""
    Lres = []
    for i in L:
        (a,b)=i
        if a == b:
            Lres.append(b)
    return len(Lres)/len(L)

    
def addBackGround(L,pbg, K):
    """Function to add background noise to a list of outcomes. With probability pbg it will add an outcome to a 
     timestep where nothing was measured
     
    Parameters:
        L : List of outcomes
        pbg: probability of getting a click due to background noise at a particular timestep
        K : Last timestep"""
    i = 0
    listestep = []
    for j in L:
        a, b = j
        listestep.append(a[0])
    while i< K:
        if i not in listestep:
            b = bernoulli.rvs(pbg)
            if b == 1:
                randbit = bernoulli.rvs(0.5)
                base = bernoulli.rvs(0.5)
                if base < 0.5:
                    L.append(([i,"plusmoins"],randbit))
                else:
                    L.append(([i,"zeroun"],randbit))
            i=i+1
        else:
            i=i+1
    L.sort(key=getTime)
    

