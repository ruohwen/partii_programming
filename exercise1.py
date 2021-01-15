import numpy as np
def mat_lin(n): #defining matrix for linear molecule with n atoms
    Mz=np.zeros((n,n)) #base matrix
    lc=Mz.shape[0] #no. of rows
    for i in range(1,lc):
        Mz[i, i-1]=-1 #atom i interacts with atom i-1
    for i in range(lc-1):
        Mz[i, i+1]=-1 #atom i interacts with atom i+1
    return(Mz)
    print(Mz)

def mat_cyc(n): #defining matrix for cyclic molecule with n atoms
    Mz=np.zeros((n,n))
    lc=Mz.shape[0]
    for i in range(1,lc):
        Mz[i, i-1]=-1
    for i in range(lc-1):
        Mz[i, i+1]=-1
    Mz[0, n-1]=-1 #atom n interacts with atom 1
    Mz[n-1, 0]=-1
    return(Mz)
    print(Mz)

def get_evals(x): #defining function to extract eigenvalues
    evals, evecs = np.linalg.eig(x)
    return(sorted(np.round(evals, 3)))
    print(sorted(np.round(evals, 3))) #print eigenvalues rounded to 3 decimal places

def degeneracies(x): #lists energy level: degeneracy
    energies=get_evals(x)
    degen = {i:energies.count(i) for i in energies}
    print(degen)

def tetrahedron():
    M=np.identity(4)-np.ones((4,4))
    return M
    print(M)

def cube():
    m = 4 #edges per face
    n = 8 #total vertices
    Mz=np.zeros((n,n))
    for i in range(0,m):
        for j in range (0,m):
            Mz[i, j]=mat_cyc(m)[i, j]
    for i in range (n-m,n,1):
        for j in range (n-m,n,1):
            Mz[i,j]=mat_cyc(m)[i-(n-m),j-(n-m)]
    for i in range(0,m):
        for j in range (n-m,n):
            Mz[i, j]=-np.identity(m)[i, j-m]
    for i in range(n-m,n):
        for j in range (0,m):
            Mz[i, j]=-np.identity(m)[i-m, j]
    return(Mz)
    print(Mz)

def dodecahedron():
    m = 5 #edges per face
    n = 20 #total vertices
    Mz=np.zeros((n,n))
    for i in range(0,m): #upper face
        for j in range (0,m):
            Mz[i, j]=mat_cyc(m)[i, j]
    for i in range (n-m,n): #lower face
        for j in range (n-m,n):
            Mz[i,j]=mat_cyc(m)[i-(n-m),j-(n-m)]
    for i in range(m,n-m): #two layers of ten points in between treated as a single polygon
        for j in range (m,n-m):
            Mz[i, j]=mat_cyc(2*m)[i-m,j-m]
    for i in range(0,5):
            Mz[i, 5+2*i]=-1
    for i in range(15,20):
            Mz[i, 2*i-24]=-1
    for i in range(0,5):
            Mz[5+2*i, i]=-1
    for i in range(15,20):
            Mz[2*i-24, i]=-1
    return(Mz)
    print(Mz)

def octahedron():
    m = 3 #edges per face
    n = 6 #total vertices
    M=np.vstack((mat_cyc(3),mat_cyc(3)))
    N=np.hstack((M,M))
    return N
    print(N)

def icosahedron():
    Mz=np.zeros((12,12))
    for i in range(1,6): #upper pentagon
        for j in range (1,6):
            Mz[i, j]=mat_cyc(5)[i-1, j-1]
    for i in range(6,11): #lower pentagon
        for j in range (6,11):
            Mz[i, j]=mat_cyc(5)[i-6, j-6]
    Mz[1:6,0]=-1 #caps
    Mz[0,1:6]=-1
    Mz[6:11,11]=-1
    Mz[11,6:11]=-1
    for i in range (1,5):
        Mz[i,i+5]=-1
        Mz[i,i+6]=-1
    for i in range (7,11):
        Mz[i,i-5]=-1
        Mz[i,i-6]=-1
    Mz[5,6]=-1
    Mz[5,10]=-1
    Mz[6,5]=-1
    Mz[6,1]=-1
    return Mz
    print(Mz)

def buckyball():
    Mz=np.zeros((60,60))
    for i in range (0,59):
        Mz[i,i+1]=-1
    Mz[0,4]=-1
    Mz[0,8]=-1
    Mz[1,11]=-1
    Mz[2,14]=-1
    Mz[3,17]=-1
    Mz[5,19]=-1
    Mz[6,21]=-1
    Mz[7,24]=-1
    Mz[9,25]=-1
    Mz[10,28]=-1
    Mz[12,29]=-1
    Mz[13,32]=-1
    Mz[15,33]=-1
    Mz[16,36]=-1
    Mz[18,37]=-1
    Mz[20,39]=-1
    Mz[22,41]=-1
    Mz[23,43]=-1
    Mz[26,44]=-1
    Mz[27,46]=-1
    Mz[30,47]=-1
    Mz[31,49]=-1
    Mz[34,50]=-1
    Mz[35,52]=-1
    Mz[38,53]=-1
    Mz[40,54]=-1
    Mz[42,56]=-1
    Mz[45,57]=-1
    Mz[48,58]=-1
    Mz[51,59]=-1
    Mz[55,59]=-1
    for i in range (0,60):
        for j in range (0,60):
            Mz[j,i]=Mz[i,j]
    return Mz
    print(Mz)