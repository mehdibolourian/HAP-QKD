import matplotlib.pyplot as plt
import numpy as np
import os

"""This scripts produces a plot of the number of successfull EPR pair transmission per second as a function of the distance between cities
 for a ballon link and a fiber link using the txt outputs from EntanglementStudy.py"""

cur_path = os.path.dirname(__file__)

path = os.path.relpath('../data', cur_path)  

NameSimu02 = os.path.join(path, "EPRFreespace.txt")         
NameTheo02 = os.path.join(path, "EPRFiber.txt")

file1 = open(NameSimu02,"r")
L = file1.read().splitlines()
key1 =[]
for i in L:
    key1.append(float(i))
file1.close()

file2 = open(NameTheo02,"r")
L = file2.read().splitlines()
key2 =[]
for i in L:
    key2.append(float(i))
file2.close()

dist_cities = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]
plt.figure(figsize=(15,9)) 
plt.semilogy(dist_cities,key1,label="With free-space links",marker = '*', markersize=15, color='m')
plt.semilogy(dist_cities,key2,label="With fiber links",marker = '+', markersize=15, color='y')

plt.xlabel('$d_{\\rm cities}$ (km) ',size=30)
plt.ylabel('Bell pairs shared per second',size=30)
plt.legend(loc='best',prop={'size':24})
plt.tick_params(axis='both', labelsize=25)
plt.savefig("EPRfiberfree.pdf", format = 'pdf')
plt.show()