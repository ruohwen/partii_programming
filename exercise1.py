import numpy as np

#defining adjacency matrix for linear molecule with n atoms
def mat_lin(n):
    Mz=np.zeros((n,n)) #zero matrix of the appropriate size
    lc=Mz.shape[0] #no. of rows
    for i in range(1,lc):
        Mz[i, i-1]=-1 #atom i interacts with atom i-1
    for i in range(lc-1):
        Mz[i, i+1]=-1 #atom i interacts with atom i+1
    return(Mz)

#defining adjacency matrix for cyclic molecule with n atoms
def mat_cyc(n):
    Mz=np.zeros((n,n))
    lc=Mz.shape[0]
    for i in range(1,lc):
        Mz[i, i-1]=-1
    for i in range(lc-1):
        Mz[i, i+1]=-1
    Mz[0, n-1]=-1 #atom n is joined to atom 1
    Mz[n-1, 0]=-1
    return(Mz)

#defining adjacency matrix of platonic solids by indexing carbons and joining them to form individual faces (similar to the cyclic polyene method). I find this more intuitive than straight-out hard coding the matrix
def drawface(M, atoms):
    atomsleft = []
    atomsright =[]
    for i in range(0,len(atoms)-1):
        atomsright.append(atoms[i+1])
    atomsright.append(atoms[0])
    atomsleft.append(atoms[len(atoms)-1])
    for i in range(1, len(atoms)):
        atomsleft.append(atoms[i-1])
    for i in range(0,len(atoms)):
        M[atoms[i]][atomsleft[i]] = -1
        M[atoms[i]][atomsright[i]] = -1
    return M

Ns = {"tetrahedron":4, "cube":8, "dodecahedron":20, "octahedron":6, "icosahedron":12} #number of vertices
facelists = {"tetrahedron": [[1,2,3],[1,3,0],[1,2,0],[2,3,0]],
             "cube": [[0,1,2,3],[4,5,6,7],[1,2,6,5],[3,2,6,7],[0,4,7,3],[0,1,5,4]],
             "dodecahedron": [[0,1,2,3,4],[5,6,7,1,0],[1,7,8,9,2],[3,2,9,10,11],[4,3,11,12,13],[0,4,13,14,5],[5,6,16,15,14],[6,7,8,17,16],[8,9,10,18,17],[10,11,12,19,18],[12,13,14,15,19]],
             "octahedron": [[0,1,4],[0,1,2],[0,2,3],[0,3,4],[5,1,4],[5,1,2],[5,2,3],[5,3,4]],
             "icosahedron": [[0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,1],[6,7,11],[7,8,11],[8,9,11],[9,10,11],[10,6,11],[2,6,7],[2,3,7],[3,7,8],[3,4,8],[4,8,9],[4,5,9],[5,9,10],[5,10,1],[1,6,10],[1,2,6]]}

def platsolfunc(selection):
    N = Ns[selection]
    facelist = facelists[selection]
    M = np.zeros((N,N))
    for face in facelist:
        atoms = face
        drawface(M, atoms)
    return M

#defining adjacency matrix of C60 using same method as platonic solids
def buckyball():
    M=np.zeros((60,60))
    facelist = [[0,1,2,3,4],
               [0,1,8,7,6,5],
               [1,2,11,10,9,8],
               [2,3,14,13,12,11],
               [3,4,17,16,15,14],
               [4,0,5,19,18,17],
               [5,6,20,39,19],
               [6,7,23,22,21,20],
               [7,8,9,24,23],
               [9,10,27,26,25,24],
               [10,11,12,28,27],
               [12,13,31,30,29,28],
               [13,14,15,32,31],
               [15,16,35,34,33,32],
               [16,17,18,36,35],
               [18,19,39,38,37,36],
               [38,39,20,21,40,54],
               [21,22,42,41,40],
               [22,23,24,25,43,42],
               [25,26,45,44,43],
               [26,27,28,29,46,45],
               [29,30,48,47,46],
               [30,31,32,33,49,48],
               [33,34,51,50,49],
               [34,35,36,37,52,51],
               [37,38,54,53,52],
               [53,54,40,41,55,59],
               [41,42,43,44,56,55],
               [44,45,46,47,57,56],
               [47,48,49,50,58,57],
               [50,51,52,53,59,58],
               [55,56,57,58,59]]
    for face in facelist:
        atoms = face
        drawface(M, atoms)
    return M

def get_evals(x): #defining function to extract eigenvalues
    evals, evecs = np.linalg.eig(x)
    return(sorted(np.round(evals, 3))) # eigenvalues rounded to 3 decimal places from most bonding to most antibonding MO

def degeneracies(x):
    energies=get_evals(x)
    degen = {i:energies.count(i) for i in energies} # creates a dictionary with key = energy value and dict[key] = degeneracy
    return(degen)

#user interface
options = ["linear polyene", "cyclic polyene", "platonic solids", "buckminsterfullerene"]
platsols = ["tetrahedron", "cube", "dodecahedron", "octahedron", "icosahedron"]
for option in options:
    print("Select " + str(options.index(option)) + " for "+ str(option))
selection = options[int(input())]
if selection == "linear polyene" or selection == "cyclic polyene":
    n = int(input("Number of atoms: "))
    if selection == "linear polyene":
        degen = degeneracies(mat_lin(n))
    elif selection == "cyclic polyene":
        degen = degeneracies(mat_cyc(n))
if selection == "platonic solids":
    for platsol in platsols:
        print("Select " + str(platsols.index(platsol)) + " for " + str(platsol))
    i = int(input())  # i for index
    selection = platsols[i]
    degen = degeneracies(platsolfunc(selection))
if selection == "buckminsterfullerene":
    degen = degeneracies(buckyball())

for key in degen:
    print("Energy level " + str(key) + " has a degeneracy of " + str(degen[key]))
print("Note: a = 0, b = -1, so the Hückel energy is the interaction energy")