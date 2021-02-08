import math
import numpy as np
import matplotlib.pyplot as plt
import time

#PROTEIN-FOLDING SIMULATION

#Extracting data from elementary reactions:
#Defining the elementary reactions: Reactants, Products, Pre-exponential factor in rate constant A, factor of concentration of urea in exp(value x [urea])
reactions=[[["D"],["I"],(2.6 * np.float_power(10, 4)), -1.68],
           [["I"],["D"],(6.0 * np.float_power(10, -2)), 0.95],
           [["I"],["N"],(7.3 * np.float_power(10, 2)), -1.72],
           [["N"],["I"],(7.5 * np.float_power(10, -4)), 1.20]]

#Initialising concentrations (arbitrary) and gradients (starting from 0 to give the get_rates() something to build upon)
Conc = {"D": 1/3, "I": 1/3, "N": 1/3}
Grad = {"D": 0, "I": 0, "N": 0}

def get_rates(concU):
    for reaction in reactions: #Going through the list of elementary reactions
        rateconstant = reaction[2] * np.exp(reaction[3] * concU) #Extract the rate constant
        for species in reaction[0]: #Going through the list of reactants
            Grad[species] += -rateconstant * Conc[species] #Add a rate of consumption of rate constant * reactant conc
        for species in reaction[1]: #Going through the list of products
            for reactant in reaction[0]:
                Grad[species] += rateconstant * Conc[reactant] #Add a rate of consumption of rate constant * REACTANT conc

concarray = [] #Create a list of lists: equilibrium concentrations of each species at different concentrations of denaturant
prevconc = {}
NewConc = {}

#Gradient descent (equilibrium): vary the concentration in time according to the rate of change until the concentration of each species plateaus
def equilibriate(concU): #takes the concentration of the denaturant as the argument
    global Grad
    if concU in range(0, 3):
        timestep = 0.000001  # forward rate constants become very large at small concU, so timestep needs to be reduced
    else:
        timestep = 0.0001 # larger step taken to allow equilibrium to be reached faster
    get_rates(concU)
    stepsize = 1
    iters = 0
    while abs(stepsize) > 0.0000001 and iters < 10000000:  # until stepsize converges to given precision or timeout
        for species in Conc:
            prevconc.update(Conc)
            NewConc[species] = Conc[species] + Grad[species] * timestep
            Conc.update(NewConc)
            stepsize = prevconc[species] - Conc[species]
        Grad = {"D": 0, "I": 0, "N": 0} #Reinitialise, recalculate gradient at new concentrations
        get_rates(concU)
        iters = iters + 1
    SumConcs = Conc["D"]+Conc["I"]+Conc["N"] #Total concentration may stray from unity
    concarray_single = []
    for species in Conc:
        concarray_single.append(Conc[species]/SumConcs) #list of equilibrium concentrations for this value of concU
    concarray.append(concarray_single) #append to the list of lists

print("Folding proteins...")

t1=time.time()
#Call equilibriate() function over [U] range
concUrange = np.arange(0,8.25,0.25)
for concU in concUrange:
    Grad = {"D": 0, "I": 0, "N": 0}
    get_rates(concU) #Initialise
    equilibriate(concU) #Self-sustaining. Calculates new conc moving forward a timestep, recalculates the rates of change at these new concentrations, and repeats until concentrations plateau.
outputconcentrationmatrix = np.vstack(concarray)
from tabulate import tabulate
print ("Equilibrium concentrations of different states of protein at varying concentrations of denaturant\n", tabulate(outputconcentrationmatrix, headers=["eq conc of D", "eq conc of I", "eq conc of N"]), "\nRows: [U]/M = ", concUrange)
t2=time.time()

#Plotting
plt.plot(concUrange, outputconcentrationmatrix[:,0], marker="x", label="D")
plt.plot(concUrange, outputconcentrationmatrix[:,1], color="red", marker="o", label="I")
plt.plot(concUrange, outputconcentrationmatrix[:,2], color="green", marker="+", label="N")
plt.title("Equilibrium Fractions of Protein States against Urea Concentration")
plt.legend()
plt.xlabel("[Urea]/M")
plt.ylabel("Fraction of Species")
plt.show()
print("Time elapsed:", round(t2-t1,2) , " seconds")

#OREGONATOR

#Extracting data, similar to above
reactions_o=[[["A","Y"],["X","P"],1.34 * np.power(10,0)],
           [["X","Y"],["P"],1.6 * np.power(10,9)],
           [["B","X"],["X","X","Z"],8 * np.power(10,3)],
           [["X","X"],["Q"],4 * np.power(10,7)],
           [["Z"],["Y"],1]]

Grad = {"A": 0,"B": 0, "P": 0, "Q": 0, "X": 0, "Y": 0, "Z": 0}
Conc = {"A": 0.06, "B": 0.06, "P": 0, "Q": 0, "X":np.float_power(10, -9.8), "Y": np.float_power(10, -6.52), "Z": np.float_power(10, -7.32)}

def get_rates_o(): #Repurposed, more general function than the one above to account for multiple reactants and multiple products, but no presence of denaturant
    for reaction_o in reactions_o:
        rateconstant = reaction_o[2]
        for species in reaction_o[0]: #for each reactant in an elementary reaction
            Grad[species] += -rateconstant * np.prod(np.array([Conc[species] for species in reaction_o[0]])) #add a rate of consumption due to this elementary reaction to the overall rate of change of this reactant. rate of consumption due to this reaction = rate constant of reaction x product of concentrations of all reactant species in this elementary reaction
        for species in reaction_o[1]: #for each product in an elementary reaction
            Grad[species] += rateconstant * np.prod(np.array([Conc[reactant_o] for reactant_o in reaction_o[0]])) #add a rate of production due to this elementary reaction. Depends on the concentration of reactant species.
#Manually coding the dependencies of each gradient on reactant concentrations actually gave a faster Oregonator

t1=time.time()
concarray = []
def oregonate(length): #Finds rate of change that is concentration-dependent, changes concentration in accordance to gradient (moving forward a timestep), recalculates gradient
    global Grad
    prevconc = {}
    NewConc = {}
    timestep = np.float_power(10, -6)
    stepsize = 1
    iters = 0
    while timestep * iters <= length: #Stop after running to a pre-defined length of time
        for species in Conc:
            prevconc.update(Conc)
            NewConc[species] = Conc[species] + Grad[species] * timestep
            Conc.update(NewConc)
            Grad = {"A": 0,"B": 0, "P": 0, "Q": 0, "X": 0, "Y": 0, "Z": 0}
            get_rates_o()
        if math.remainder(iters, 1000) == 0: #Plot one dot every 1000
            concarray_pointintime = [iters*timestep, Conc["X"], Conc["Y"], Conc["Z"]]
            concarray.append(concarray_pointintime)
        iters = iters + 1

print("Oregonator is oregonating...")
get_rates_o()
oregonate(0.2) #Define the length of the time axis here
output = np.vstack(concarray)
np.savetxt("oregonator.csv", output, delimiter=",")
print(output)
t2 = time.time()
print("Time elapsed: ", round((t2-t1)/60, 2), " minutes")

timeseries = output[:,0]

#Plotting
plt.plot(timeseries, output[:,1], label="X")
plt.plot(timeseries, output[:,2], color="red", label="Y")
plt.plot(timeseries, output[:,3], color="green", label="Z")
plt.yscale("log")
plt.title("Concentrations of Species X, Y, and Z in the Oregonator over time")
plt.legend()
plt.xlabel("t/s")
plt.ylabel("Concentration/M (logarithmic)")
plt.show()

'''
Spatially-resolved simulation of the Oregonator:
- Define a class called Cell which stores the concentrations of all species and contains methods to obtain rates of change of each species, take a step forward in time (vary concentration), and to interact with its neighbouring cells by diffusion
- In 1D, split a pot up into 11 cells indexed from -5 to 5. Initialise Cell(0) with the starting concentrations of reactants (all other cells, all concentrations 0)
- Defining diffusion: At time t, concentration of species x in any cell i containing it is reduced to a third of its original value. The adjacent cells (i+1 and i-1) gain a third each of the concentration in cell i. 
- Repeat for all cells containing species x, except for cells -5 and 5 (next to the walls of the pot), where the concentration is reduced by half and propagated to cells -4 and 4 only.
- Need to define a diffusion rate (time values at which the diffusion method is called). Can be defined for individual species.
- Each cell running an Oregonator internally, constantly updating the concentration independently of diffusion.
- In 2D (square cells), Cell(i,j) interacts with Cells(i+1,j),(i-1,j),(i,j-1),(i,j+1). Concentration divided over five cells for each diffusion step.
- Each cell outputs concentrations of X, Y, and Z in time. Can assign a colour to each (RGB) to visualise. Plotting intensity against position at a given t value gives a snapshot.
- Need to have specific initial concentrations for any colour-changing patterns to be observed at all.
'''