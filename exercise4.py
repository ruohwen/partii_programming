import numpy as np
import math
from random import randint
import time
import sys

class vec3d:  #Defining 3D vectors as an object type
    def __init__(self, lst = (0,0,0)):
        self.coords = list(lst)

    def magnitude(self):
        return math.sqrt(self.coords[0]*self.coords[0] + self.coords[1]*self.coords[1] + self.coords[2]*self.coords[2])

    def __sub__(self, other): #defining vector subtraction as new x coordinate = difference between x coordinate of vector 1 and x coordinate of vector 2
        ret = vec3d()
        for i in range (3):
            ret.coords[i] = self.coords[i] - other.coords[i]
        return ret

#User input: number of particles; initialise positions
n = int(input("Number of particles in the system: "))
positionupperlimit = (n*10)/2 #Defining a reasonable range of starting positions so that energies converge, which depends on the number of particles in the system
r_dummy = []
r = []
for x in range (0,n):
    r_dummy.append(vec3d((randint(0,positionupperlimit),randint(0,positionupperlimit),randint(0,positionupperlimit))).coords) #Random positions within reasonable range defined above, to initialise.
for x in range(0,n):
    r.append([float(item) for item in r_dummy[x]]) #Coordinates in vector class are a list of integers. Copy out into r (from r_dummy) and change type into float to allow positions in steps of 0.1.
r = np.vstack(r)
r = r/10 #Allow positions to be initialised in steps of 0.1 (the chances of overlapping particles were too high when positions were defined in steps of 1)
print("Randomly-generated position vectors for ", n, " particles")
print(r)

#Define potential function and energy unit. Allows the user to specify a general potential in addition to the pre-defined LJ and Morse potentials.
PEselection = input("Potential curve\nSelect 1 for Lennard-Jones potential\nSelect 2 for Morse potential\nSelect 3 to specify your own function\n")
if PEselection == "1":
    def u(R):
        return 4 * (1/np.power(R, 12)-1/np.power(R, 6))
    energyunit = "ε"
elif PEselection == "2":
    R_e = float(input("Input value of (r_e)/σ: "))
    def u(R):
        return np.power((1-np.exp(-(R-R_e))), 2)
    energyunit = "D_e"
elif PEselection == "3":
    inputfunction = input("Specify your pairwise interaction potential as a function of R (unit σ): ")
    def u(R):
        return eval(inputfunction)
    energyunit = input("Unit of energy: ")
else:
    print("Invalid selection. Please try again.")
    sys.exit(0)

#Function to calculate the total potential energy of the configuration (defined by coordinates in r and potential function u(R)), by summing all pairwise interaction energies
def U(r):
    U = 0
    for i in range (0,n):
        for j in range (i+1, n):
            R = (vec3d((tuple(r[j])))-vec3d((tuple(r[i])))).magnitude()
            U += u(R)
    return U

U(r)
print("The energy of this initial system is", U(r), energyunit)

#Computing energy gradient using first principles

h = 0.001 #Step size, found this to be a reasonable value by trial and error. Should be infinitismally small
r_forward = r.copy()
r_backward = r.copy()

def d_u(r, i, j): #i is particle index (runs from 0 to n-1), j is coordinate index i.e. for (x, y, z), j=(0, 1, 2)
    r_forward[i][j] = r[i][j] + h #takes step forward in j direction. New position vectors.
    r_backward[i][j] = r[i][j] - h #takes step backward in j direction. New position vectors.
    U_forward = 0
    U_backward = 0
    for p in range (0,n):
        for q in range (p+1, n):
            R_forward = (vec3d((tuple(r_forward[q])))-vec3d((tuple(r_forward[p])))).magnitude() #Gives the distance between particle p and particle q in the new coordinates r_forward, i.e. when particle i has moved 0.001 in the j direction
            U_forward += u(R_forward) #Calculates the PE of the system in the limit of a small perturbation in the position of particle i in the j direction
    for p in range (0,n):
        for q in range (p+1, n):
            R_backward = (vec3d((tuple(r_backward[q])))-vec3d((tuple(r_backward[p])))).magnitude() #Gives the distance between particle p and particle q in the new coordinates r_backwardi.e. when particle i has moved 0.001 in the -j direction
            U_backward += u(R_backward)#Calculates the PE of the system in the limit of a small perturbation in the position of particle i in the -j direction
    r_forward[i][j] = r[i][j]  # return to original positions
    r_backward[i][j] = r[i][j]
    return (U_forward - U_backward) / (2*h) #symmetrical energy gradient when particle i moves in j direction

#Function to move particle i a small step down the energy slope along j
def step(r, i, j):
    if d_u(r, i, j) < 0:
        r[i][j] = r[i][j] + h
    elif d_u(r, i, j) > 0: #if +j increases U, move in -j
        r[i][j] = r[i][j] - h
    else: #gradient = 0, already occupying minimum along j. Unlikely.
        r[i][j] = r[i][j]
    return r

#Function for gradient descent until minimum, varying the coordinate of a single particle i in a single direction j. All other coordinates and particles unchanged.
def descend(r, i, j):
    iters = 0
    d_u(r, i, j)
    while iters < 3000 and abs(d_u(r,i,j)) >= 0.1:
        d_u(r,i,j) #Samples gradient
        step(r,i,j) #Take a step in the direction which lowers energy
        iters = iters + 1 #Repeat until gradient flattens out sufficiently (smaller than 0.1) to the energy minimum, or aborts if energy minimum is still not located after moving particle i more than 3σ in j direction (fine because will repeat and take global minimum)
    print("Energy =",U(r), energyunit)

#"Run": Execute gradient descent for each particle in each direction, in turn. Repeat until energy converges (a local minimum is reached in all i and j coordinates).
Optimisation=[]
def run():
    for i in range(n):
        for j in range(3):
            descend(r, i, j)
    Optimisation.append(U(r))

list_energies = []
list_configurations = []
list_distances = []
list_iters = []

def singlerun(): #execute run() to find a local minimum, append the resulting energy and configuration to lists
    distances = []
    run()
    run()
    iters = 2
    while abs(Optimisation[-1]-Optimisation[-2])>0.00001: #Run gradient descent (all particles, all directions) until energies converge to a precision as given
        run()
        iters = iters+1
    list_energies.append(round(Optimisation[-1],2))
    list_configurations.append(r)
    for i in range(0, n):
        for j in range(i + 1, n):
            distances.append([int(i), int(j), ((vec3d((tuple(r[j]))) - vec3d((tuple(r[i])))).magnitude())])
    list_distances.append(np.vstack(distances))
    list_iters.append(iters)

#Run whole program a user-defined number of times in an attempt to find the global minimum
noofattempts = int(input("Please input the number of runs to attempt:\n(More runs are recommended for larger systems and especially when using the LJ potential to find the global minimum, but please expect a wait of ~1 min per run.\nRecommended 5 runs for LJ potential)\n"))
t1 = time.time()
print("Optimising...")
for runnumber in range(0,noofattempts):
    print("Run "+str(runnumber+1)+" start")
    singlerun()
    print("Run " + str(runnumber+1) + " end")
    if runnumber+1 != noofattempts:
        print("Reinitialising")
        r_dummy = []
        for x in range(0, n):
            r_dummy.append(vec3d((randint(0, positionupperlimit), randint(0, positionupperlimit), randint(0,positionupperlimit))).coords)  # Random positions within reasonable range defined above, to initialise.
        r = []
        for x in range(0, n):
            r.append([float(item) for item in r_dummy[x]])
        r = np.vstack(r)
        r = r / 10
        print("Randomly-generated position vectors for ", n, " particles")
        print(r)

#Identifying configutarion corresponding to global energy minimum by parsing through the list generated
lowestenergy = min(list_energies)
lowestenergy_runnumber = list_energies.index(lowestenergy)
optimalconfig = list_configurations[lowestenergy_runnumber]
optimaldistances = list_distances[lowestenergy_runnumber]
optimaliters = list_iters[lowestenergy_runnumber]

t2 = time.time()

#OUTPUT
#Potential energy (global minimum), new coordinates, inter-particle distances
print("Energies obtained in " +str(noofattempts)+" runs:"+ str(list_energies))
print("Lowest energy of", lowestenergy, energyunit, "was obtained on run", lowestenergy_runnumber+1, "after", optimaliters, "iterations.")
from tabulate import tabulate
print("The positions of the particles at this lowest energy conformation is\n"+tabulate(optimalconfig, headers=["x", "y", "z"]))
print ("Inter-particle distances\n"+tabulate(optimaldistances, headers=["Particle i", "Particle j", "Distance between i and j / σ"]))

#XYZ file
f = open("output.xyz", "w")
f.write(str(n))
f.write("\nCoordinates of particles in minimum energy configuration\n")
for x in range(n):
    line = str("H\t" + str(round(optimalconfig[x][0],5)) + "\t"+ str(round(optimalconfig[x][1],5)) + "\t" + str(round(optimalconfig[x][2],5)) + "\t\n")
    f.write(line)
f.close()

print("Time elapsed:")
if t2-t1 >= 60:
    print(round((t2-t1)/60, 2), "minutes")
else:
    print(round(t2-t1, 2), "seconds")