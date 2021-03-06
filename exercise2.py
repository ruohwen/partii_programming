import os
import numpy as np

#Extracting energy values for each bond length and angle
d1={}
d2={}
d3={}
directory = input("Enter directory path: ")
print("Parsing and plotting the data may take some time, please wait...")
files=os.listdir(directory) #parse directory

for file in files:
    f=open(os.path.join(directory,file),'r')
    lines=f.readlines()
    tls=list(file) #split file name into individual elements
#In hindsight, it's probably better to take the bond lengths and angles from inside the file rather than the file name, in case the file names get changed. But it seemed like the easiest option at the time.
    bondlength=float(str(tls[5])+str(tls[6])+str(tls[7])+str(tls[8]))
    angle=str(tls[-8])+str(tls[-7])+str(tls[-6])+str(tls[-5])
    if str(tls[-9]) == "1": #for three digit bond angles (>=100.0)
        angle=str(tls[-9])+str(angle)
    angle=float(angle)
    for i,line in enumerate(lines): #parse for energy
        if "SCF" in line:
            energyline=lines[i-5]
    words = str.split(energyline)
    energy = float(words[4])
    d1[file] = bondlength #building dictionaries for bond length, angle, and energy respectively, each with key = filename
    d2[file] = angle
    d3[file] = energy
X=[]
Y=[]
Z=[]
for file in files:
    X.append(d1.get(file)) #creating 1D arrays for bond length, angle, and energy (call from dictionary)
    Y.append(d2.get(file))
    Z.append(d3.get(file))

#Plotting potential energy surface
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(np.array(X), np.array(Y), np.array(Z), cmap=cm.magma)
ax.set_xlabel('r/Å')
ax.set_ylabel('θ/°')
ax.set_zlabel('Energy/Hartrees')
ax.set_title("Potential energy surface for molecule in directory " + str(directory))

plt.show()

#Locating PE minimum
energymin = min(Z)
zindex = Z.index(energymin) #Finding the lowest energy value in the energy list, take its index, and extract the corresponding bond length and angle
bondlengthformin = X[zindex]
angleformin = Y[zindex]
print("Energy minimum (r/Å, θ/°, Energy/Hartrees) is observed at ("+str(bondlengthformin)+", "+str(angleformin)+", "+str(energymin)+")")

from scipy.optimize import curve_fit
matrix=np.ndarray.transpose(np.vstack((X,Y,Z)))

#Fitting curves around energy minimum
print("Defining range around energy minimum to fit to curve:")
bondlengthrange=input("±Range of bond length (recommended: 0.1)/Angstrom: ")
anglerange=input("±Range of angles (recommended:9)/degrees: ")
#The area taken is a square, although it might be preferable to take a circle. Molecules with bond length and angle close to that at the energy minimum preferred (over purely taking lowest 100 or so energies) in case there are several potential minima - the graph seems to show that this is not the case (but we cannot be sure)
neareqbl=[] #create lists of near equilibrium bond lengths, angles, and energies
neareqang=[]
neareqenergy=[]
for i in range(matrix.shape[0]):
    if matrix[i,1] in range (int(angleformin) - int(anglerange), int(angleformin) + int(anglerange)):
        if (100*matrix[i, 0]) in range(int(100*(bondlengthformin - float(bondlengthrange))), int(100*(bondlengthformin + float(bondlengthrange)))):
            neareqbl.append(matrix[i,0])
            neareqang.append(matrix[i,1])
            neareqenergy.append(matrix[i,2])
displacement=np.array(neareqbl) - bondlengthformin
angdisplacement=np.array(neareqang) - angleformin
hartreetojoule = 4.3597*pow(10, -18)
neareqenergy = np.array(neareqenergy) * hartreetojoule

def func(X,kr,ktheta):
    X = np.array([displacement,angdisplacement])
    return energymin*hartreetojoule+(1/2)*kr*np.square(displacement*pow(10,-10))+(1/2)*ktheta*np.square(angdisplacement*(2*np.pi/360))

popt, pcov = curve_fit(func, [displacement,angdisplacement], neareqenergy)
kr=popt[0]
ktheta=popt[1]
print("Taking "+str(len(neareqenergy))+" points around the minimum, kr is "+str(kr) +"Nm^-1 and ktheta is "+ str(ktheta)+ "Jrad^-2")

#Calculating vibrational frequencies
mu=1.66*pow(10,-27)
v1= (1/(2*np.pi)) * np.sqrt(kr/(2*mu))
v2= (1/(2*np.pi)) * np.sqrt(ktheta/(np.square(bondlengthformin*pow(10,-10))*0.5*mu))
Hztopercm=3.3356*pow(10,-11)
v1percm = v1*Hztopercm
v2percm = v2*Hztopercm
print("The vibrational frequencies for the stretching and bending modes, v1 and v2, are "+str(v1percm)+"cm^-1 and "+str(v2percm)+"cm^-1 respectively.")