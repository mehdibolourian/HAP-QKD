import matplotlib.pyplot as plt
import numpy as np
import os

"""This scripts plots the transmissivity of the vertical Balloon-To-Ground downlink channel using the txt outputs from HeightStudy.py"""
cur_path = os.path.dirname(__file__)

path = os.path.relpath('../data', cur_path)  

NameSimu02 = os.path.join(path, "ZenithBalloonSimu01.txt")         
NameTheo02 = os.path.join(path, "ZenithBalloonTheo01.txt")
NameSimu03 = os.path.join(path, "ZenithBalloonSimu02.txt")
NameTheo03 = os.path.join(path, "ZenithBalloonTheo02.txt")
NameSimu04 = os.path.join(path, "ZenithBalloonSimu04.txt")
NameTheo04 = os.path.join(path, "ZenithBalloonTheo04.txt")
NameSimu05 = os.path.join(path, "ZenithBalloonSimu06.txt")
NameTheo05 = os.path.join(path, "ZenithBalloonTheo06.txt")
NameSimu06 = os.path.join(path, "ZenithBalloonSimu08.txt")
NameTheo06 = os.path.join(path, "ZenithBalloonTheo08.txt")

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

file3 = open(NameSimu03,"r")
L = file3.read().splitlines()
key3 =[]
for i in L:
    key3.append(float(i))
file3.close()

file4 = open(NameTheo03,"r")
L = file4.read().splitlines()
key4 =[]
for i in L:
    key4.append(float(i))
file4.close()

file5 = open(NameSimu04,"r")
L = file5.read().splitlines()
key5 =[]
for i in L:
    key5.append(float(i))
file5.close()

file6 = open(NameTheo04,"r")
L = file6.read().splitlines()
key6 =[]
for i in L:
    key6.append(float(i))
file6.close()

file7 = open(NameSimu05,"r")
L = file7.read().splitlines()
key7 =[]
for i in L:
    key7.append(float(i))
file7.close()

file8 = open(NameTheo05,"r")
L = file8.read().splitlines()
key8 =[]
for i in L:
    key8.append(float(i))
file8.close()

file9 = open(NameSimu06,"r")
L = file9.read().splitlines()
key9 =[]
for i in L:
    key9.append(float(i))
file9.close()

file10 = open(NameTheo06,"r")
L = file10.read().splitlines()
key10 =[]
for i in L:
    key10.append(float(i))
file10.close()

dist = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,77,80]
y_error = 0.002

plt.figure(figsize=(15,9)) 
plt.plot(dist,key1,linestyle = '', marker = '*', markersize=10, color='m')
plt.plot(dist,key2,label="$H=18$ km", color='m')

plt.plot(dist,key3,linestyle = '', marker = '*', markersize=10, color='g')
plt.plot(dist,key4,label="$H=23$ km", color='g')

plt.plot(dist,key5,linestyle = '', marker = '*', markersize=10,color='r')
plt.plot(dist,key6,label="$H=28$ km",color='r')

plt.plot(dist,key7,linestyle = '', marker = '*', markersize=10,color='b')
plt.plot(dist,key8,label="$H=33$ km",color='b')

plt.plot(dist,key9,linestyle = '', marker = '*', markersize=10,color='y')
plt.plot(dist,key10,label="$H=38$ km",color='y')

plt.xlabel('Zenith angle (degree) ',size=30)
plt.ylabel('Mean channel efficiency',size=30)
plt.legend(loc='upper right',prop={'size':24})
plt.tick_params(axis='both', labelsize=25)

plt.savefig("ZenithHeight.pdf", format = 'pdf')
plt.show()