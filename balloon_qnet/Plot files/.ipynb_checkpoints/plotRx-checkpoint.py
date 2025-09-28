import matplotlib.pyplot as plt
import numpy as np
import os
"""This scripts plots the transmissivity of the vertical Balloon-To-Ground downlink channel using the txt outputs from RxStudy.py"""

cur_path = os.path.dirname(__file__)

path = os.path.relpath('../data', cur_path)  

NameSimu02 = os.path.join(path, "HeightballoonSimu00.txt")         
NameTheo02 = os.path.join(path, "HeightballoonTheo00.txt")
NameSimu03 = os.path.join(path, "HeightballoonSimu01.txt")
NameTheo03 = os.path.join(path, "HeightballoonTheo01.txt")
NameSimu04 = os.path.join(path, "HeightballoonSimu02.txt")
NameTheo04 = os.path.join(path, "HeightballoonTheo02.txt")
NameSimu05 = os.path.join(path, "HeightballoonSimu04.txt")
NameTheo05 = os.path.join(path, "HeightballoonTheo04.txt")
NameSimu06 = os.path.join(path, "HeightballoonSimu06.txt")
NameTheo06 = os.path.join(path, "HeightballoonTheo06.txt")

file0 = open(NameSimu02,"r")
L = file0.read().splitlines()
key0 =[]
for i in L:
    key0.append(float(i))
file0.close()

file01 = open(NameTheo02,"r")
L = file01.read().splitlines()
key01 =[]
for i in L:
    key01.append(float(i))
file01.close()

file1 = open(NameSimu03,"r")
L = file1.read().splitlines()
key1 =[]
for i in L:
    key1.append(float(i))
file1.close()

file2 = open(NameTheo03,"r")
L = file2.read().splitlines()
key2 =[]
for i in L:
    key2.append(float(i))
file2.close()

file3 = open(NameSimu04,"r")
L = file3.read().splitlines()
key3 =[]
for i in L:
    key3.append(float(i))
file3.close()

file4 = open(NameTheo04,"r")
L = file4.read().splitlines()
key4 =[]
for i in L:
    key4.append(float(i))
file4.close()

file5 = open(NameSimu05,"r")
L = file5.read().splitlines()
key5 =[]
for i in L:
    key5.append(float(i))
file5.close()

file6 = open(NameTheo05,"r")
L = file6.read().splitlines()
key6 =[]
for i in L:
    key6.append(float(i))
file6.close()

file7 = open(NameSimu06,"r")
L = file7.read().splitlines()
key7 =[]
for i in L:
    key7.append(float(i))
file7.close()

file8 = open(NameTheo06,"r")
L = file8.read().splitlines()
key8 =[]
for i in L:
    key8.append(float(i))
file8.close()
dist = range(18,38)
y_error = 0.002

plt.figure(figsize=(12,10))
plt.plot(dist,key0, linestyle = '', marker = '*', markersize=10, color='y')
plt.plot(dist,key01,label="$D_{\\mathrm{Rx}}=20$ cm", color='y')

plt.plot(dist,key1, linestyle = '', marker = '*', markersize=10, color='m')
plt.plot(dist,key2,label="$D_{\\mathrm{Rx}}=30$ cm", color='m')

plt.plot(dist,key3, linestyle = '', marker = '*', markersize=10, color='g')
plt.plot(dist,key4,label="$D_{\\mathrm{Rx}}=40$ cm", color='g')

plt.plot(dist,key5, linestyle = '', marker = '*', markersize=10,color='r')
plt.plot(dist,key6,label="$D_{\\mathrm{Rx}}=50$ cm",color='r')

plt.plot(dist,key7, linestyle = '', marker = '*', markersize=10,color='b')
plt.plot(dist,key8,label="$D_{\\mathrm{Rx}}=60$ cm",color='b')

# plt.errorbar(dist, key1,yerr = y_error,fmt ='*',ecolor='m')
# plt.errorbar(dist, key3,yerr = y_error,fmt ='*',ecolor='g')
# plt.errorbar(dist, key5,yerr = y_error,fmt ='*',ecolor='r')
# plt.errorbar(dist, key7,yerr = y_error,fmt ='*',ecolor='b')

plt.xlabel('Height (km) ',size=30)
plt.ylabel('Mean channel efficiency',size=30)
plt.legend(loc='best',prop={'size':28})
plt.tick_params(axis='both', labelsize=28)

plt.savefig("HeightStudy.pdf", format = 'pdf')
plt.show()