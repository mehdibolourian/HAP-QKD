import matplotlib.pyplot as plt
import numpy as np
import os

"""This scripts produces a plot of the transmissivity of the vertical Balloon-To-Ground downlink channel for different initial beam waists
 using the txt outputs from W0Study.py"""

cur_path = os.path.dirname(__file__)

path = os.path.relpath('../data', cur_path)  

NameSimu01 = os.path.join(path, "W0Simu01.txt")         
NameTheo01 = os.path.join(path, "W0Theo01.txt")
NameSimu02 = os.path.join(path, "W0Simu02.txt")
NameTheo02 = os.path.join(path, "W0Theo02.txt")
NameSimu03 = os.path.join(path, "W0Simu03.txt")
NameTheo03 = os.path.join(path, "W0Theo03.txt")
NameSimu04 = os.path.join(path, "W0Simu04.txt")
NameTheo04 = os.path.join(path, "W0Theo04.txt")

file1 = open(NameSimu01, "r")
L = file1.read().splitlines()
key1 =[]
for i in L:
    key1.append(float(i))
file1.close()

file2 = open(NameTheo01,"r")
L = file2.read().splitlines()
key2 =[]
for i in L:
    key2.append(float(i))
file2.close()

file3 = open(NameSimu02,"r")
L = file3.read().splitlines()
key3 =[]
for i in L:
    key3.append(float(i))
file3.close()

file4 = open(NameTheo02,"r")
L = file4.read().splitlines()
key4 =[]
for i in L:
    key4.append(float(i))
file4.close()

file5 = open(NameSimu03,"r")
L = file5.read().splitlines()
key5 =[]
for i in L:
    key5.append(float(i))
file5.close()

file6 = open(NameTheo03,"r")
L = file6.read().splitlines()
key6 =[]
for i in L:
    key6.append(float(i))
file6.close()

file7 = open(NameSimu04,"r")
L = file7.read().splitlines()
key7 =[]
for i in L:
    key7.append(float(i))
file7.close()

file8 = open(NameTheo04,"r")
L = file8.read().splitlines()
key8 =[]
for i in L:
    key8.append(float(i))
file8.close()

dist = range(18,38)

plt.figure(figsize=(15,9)) 
plt.plot(dist,key1, linestyle = '', marker = '*', markersize=10, color='m',alpha=0.5)
plt.plot(dist,key2,label="$W_0 = 5$ cm", color='m')

plt.plot(dist,key3, linestyle = '', marker = '*', markersize=10, color='g')
plt.plot(dist,key4,label="$W_0 = 10$ cm", color='g')

plt.plot(dist,key5, linestyle = '', marker = '*', markersize=10,color='r')
plt.plot(dist,key6,label="$W_0 = 15$ cm",color='r')

plt.plot(dist,key7, linestyle = '', marker = '*', markersize=10,color='b')
plt.plot(dist,key8,label="$W_0 = 20$ cm",color='b')

plt.xlabel('Height (km) ',size=30)
plt.ylabel('Mean channel efficiency',size=30)
plt.legend(loc='best',prop={'size':24})
plt.tick_params(axis='both', labelsize=25)
plt.savefig("W0Trans.pdf", format = 'pdf')
plt.show()