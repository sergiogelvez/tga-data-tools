import pandas as pd
import scipy.integrate as integrate
import scipy.special as special
from math import exp, inf, log
import numpy as np
import matplotlib.pyplot as pp
#import matplotlib.axes as ax

# acá se cambia el camino
inputFile = open("PMMA.csv")
outputFile = open("PMMA-salida.csv",'w')

T = np.array([])
E = np.array([])


for line in inputFile:
    #print " ".join(line.split())
    #print line
    l = line.split(";")
    #print linecita
    T = np.append(T, float(l[1]))
    E = np.append(E, float(l[4]))
    #lC.append(float(l[2]))
    
    
inputFile.close()
#file2.close()

#E = 320000
R = 8.314
#T = np.array([773.15, 900])
#T = np.linspace(873.15, 298.15, 10)
#E = np.array([270.977,270.935,270.893,270.85,270.808,270.766,270.724,270.682,270.64,270.598,270.556,270.514,270.472,270.43,270.388,270.346,270.304,270.262,270.22,270.178,270.136,270.094,270.052])
#T = np.array([512.85,512.847,512.844,512.84,512.837,512.834,512.83,512.827,512.824,512.82,512.817,512.814,512.81,512.807,512.804,512.8,512.797,512.794,512.79,512.787,512.784,512.78,512.777])
#E = np.array([266.332,266.291,266.249,266.208,266.167,266.126,266.085,266.043,266.002,265.961])
E = E * 1000
T = T + 273.15
#T = np.array([512.48,512.477,512.474,512.47,512.467,512.464,512.46,512.457,512.454,512.45])
print("Valor de T: "+str(T))
beta = 2
Z = np.array([])
for ei, ti in zip(E,T):
    zi = ei/(R*ti)
    Z = np.append(Z, zi)

print("Valor de E: " + str(E))
# 
print("Zs: "+str(Z))

# crear P para guardar los p para cada T
P = np.array([])

#variable para controlar el numero de pasos que se han ejecutado, para depurar
i = 0

for zi in Z:
    print("paso numero " + str(i))
    try:
        result = integrate.quad(lambda x: 1/(x*exp(x)), zi, +inf)
    except OverflowError:
        print(zi)
        #result[0] = float('inf')
        #result[1] = float('inf')
    #print(result)
    #result2 
    pi = 1/(exp(zi)*zi) - result[0]
    print("Valor de p(xi) = " + str(pi))
    P = np.append(P, pi)
    i = i + 1

print("Hizo todos los calculos de p(x)")
#print(P)

ltf = np.array([])
tf = np.array([])
for ti, pi, ei in zip(T, P, E):
    ltfi = ei/(R*ti) + log((ei/(beta*R))*pi)
    print(ltfi)
    print(1000/ti)
    ltf = np.append(ltf, ltfi)
    tf = np.append(tf, exp(ltfi))

print("ltf: ")
print(ltf)
print("tf: ")
print(tf)

print("tamaño de tf: " + str(tf.size))
print("tamaño de ltf: " + str(ltf.size))
print("tamaño de T: " + str(T.size))
print("tamaño de E: " + str(E.size))

# escribir un archivo con columnas T, 1000/T, ltf, tf
for ti, ltfi, tfi in zip(T, ltf, tf):
    outputFile.write(str(ti) + "; " + str(1000/ti) + "; " + str(ltfi) + "; " + str(tfi) + "\n")

outputFile.close()

#pp.plot(1000/T, ltf)
pp.plot(T, E)
#pp.plot(1000/T, tf)
pp.gca().invert_xaxis()
pp.show()
