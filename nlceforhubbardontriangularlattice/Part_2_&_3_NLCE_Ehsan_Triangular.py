# This script was first developed by Jyoti Rani in 2021 and later modified, debugged, and further documented by Ehsan Khatami in 2022.


import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from itertools import permutations 


# 
# Input: NLCE_1_#.txt stores the graphs for Order #
# 
# Outputs: NLCE_1_#_data.npy stores the Ising model data for Order # as a binary file 
# 

# In this part, we read in the bond information for every 
# cluster and solve the Ising model for them.

# Maximum NLCE order
Nmax = 8

# Number of temperatures
NT = 100

# Choosing a logarithmic temperature grid
Temp = np.logspace(-1,2,NT)

# Holding "lattice constants" (multiplicities) for all clusters
# in all orders. Probably good for up to order 9. Need to check 
# what is needed exactly for the size of the second dimension.
multi = np.zeros([Nmax+1, 3500], int)

# Holding the two site numbers and optionally the bond type 
# for every bond.
site1 = np.zeros(40,int)
site2 = np.zeros(40,int)
typ   = np.zeros(40,int)

# Exchange constant and the magnetic field
J = np.array([1, 0])
Bfield = 0.0

# This function helps determine if a bit represents 
# spin-up or spin-down
def isKthBitSet(n, k):
    if n & (1 << k):
        return 1
    else:
        return -1

# Counter for isomorphically distinct clusters
topN = 0

# Loop over the NLCE order
for N in range(2, Nmax+1):
    
  # Initializing arrays and openning files
  Estore = np.zeros([1500,NT])
  Mstore = np.zeros([1500,NT])
  Esqstore = np.zeros([1500,NT])
  Msqstore = np.zeros([1500,NT])

  # Change 1 to 2 in the line below and in the next cell 
  # to include n.n.n. bonds
  fbase = "NLCE_1_" + str(N)
  fname =  fbase + ".txt"
  file = open(fname, 'r')
  fnameE = fbase + "_dataE.npy"
  fnameM = fbase + "_dataM.npy"
  fnameEsq = fbase + "_data_Esq.npy"
  fnameMsq = fbase + "_data_Msq.npy"
  fileE = open(fnameE, 'wb')
  fileM = open(fnameM, 'wb')
  fileEsq = open(fnameEsq, 'wb')
  fileMsq = open(fnameMsq, 'wb')
 
  # Skips line with order number 
  next(file)

  # Going through the file for each order line by line
  # and reading in the cluster information
  topN = int(file.readline())
  print("ORDER", N)

  EOF = False
  while EOF == False:

    line = file.readline().split()

    # Get the number of bonds from this line
    nB = int(line[0])

    # Read in the bond information
    for b in range(nB):
      line = file.readline().split()
      site1[b] = int(eval(line[0]))
      site2[b] = int(eval(line[1]))
      typ[b]   = int(eval(line[2])) -1

# Doing the Ising model and storing four quantieis:
#
# 1- Average energy, <H>
# 2- Average magnetization, <M>
# 3- Average energy squared, <H^2>, and
# 4- Average magnetization squared, <M^2>
#
# The last two quantities will be used to compute the heat
# capacity and magnetic susceptibility, respectively. See
# Eq. (7) of https://arxiv.org/pdf/1204.1556.pdf for the 
# heat capacity of the Hubbard model in terms of 
# correlation functions.
#---------------------------------------------------------
    numerator_E   = np.zeros(NT)
    numerator_M   = np.zeros(NT)
    numerator_Esq = np.zeros(NT)
    numerator_Msq = np.zeros(NT)
    denominator   = np.zeros(NT)

    for n in range(2**N):
        
        # Magnetization and energy for each configuration
        M_n = 0
        E_n = 0
        
        # We use bit operations to determine spin directions
        for k in range(N):
            M_n += isKthBitSet(n,k)
        for b in range(nB):
            E_n += J[typ[b]]*isKthBitSet(n,site1[b])*isKthBitSet(n,site2[b])
            
        # Optionally add contribution from the magentic field
        E_n -= M_n*Bfield

        # Thermal sums are done here
        for iT in range(NT):
            P_n = np.exp(-E_n/Temp[iT])
            numerator_M[iT]   += M_n*P_n
            numerator_E[iT]   += E_n*P_n
            numerator_Esq[iT] += E_n**2*P_n
            numerator_Msq[iT] += M_n**2*P_n
            denominator[iT]   += P_n

    Estore[topN,:] = numerator_E/denominator
    Mstore[topN,:] = numerator_M/denominator
    
    # It is important to do the following subtractions at the 
    # individual cluster level to make the quantities extensive
    # for the NLCE.
    Esqstore[topN,:] = numerator_Esq/denominator - Estore[topN,:]**2
    Msqstore[topN,:] = numerator_Msq/denominator - Mstore[topN,:]**2
#---------------------------------------------------------

# Here, we take the opportunity to read in the "lattice constants" 
# (multiplicities) for each topological cluster if we have reached 
# the end of the file.

    # skipping the subgraph information for now
    sc_Num = int(file.readline())
    for s in range(sc_Num):
      next(file)
    next(file)
    
    # Checking if we have reached the end of the file
    if "Topo. # : Multp." in file.readline():
      #next(file)
      for i in range(topN+1):
        multi[N][i] = file.readline().split()[1]
      EOF = True
      file.close()
      break
    else:
      # If not yet done with clusters 
      topN += 1

  # Saving the properties to files
  np.save(fileE, Estore)
  np.save(fileM, Mstore)
  np.save(fileEsq, Esqstore)
  np.save(fileMsq, Msqstore)

  fileE.close()
  fileM.close()
  fileEsq.close()
  fileMsq.close()



# In this part, we use the ED results, the subgraph information
# and the lattice constants and construct the NLCE sums.

# This array is going to hold the contribution of clusters 
# from all orders
weights = np.zeros([4, Nmax+1, 3500, NT])

# Cluster counter; this is the same as topN from the previous part
c = 0

def NLCE():
  # O is going to hold NLCE's partial sums (sums involving clusters
  # in a particular order order) for our four properties. These are
  # Sn in Eq. (3) of https://arxiv.org/pdf/0706.3254.pdf
  O = np.zeros([4,Nmax+1,NT])

  # Hard coding the contributions from the single site by hand
  singleE = np.zeros(NT)
  singleM = np.tanh(Bfield/Temp) 
  singleEsq = np.zeros(NT)
  singleMsq = np.ones(NT)

  # Loop over NLCE orders
  for N in range(2, Nmax+1):
    
    # Opening the property binary files and loading the information
    # on to arrays
    fbase = "NLCE_1_" + str(N)
    fname =  fbase + ".txt"
    
    file = open(fname, 'r')

    fnameE = fbase + "_dataE.npy"
    fnameM = fbase + "_dataM.npy"
    fnameEsq = fbase + "_data_Esq.npy"
    fnameMsq = fbase + "_data_Msq.npy"
    fileE = open(fnameE, 'rb')
    fileM = open(fnameM, 'rb')
    fileEsq = open(fnameEsq, 'rb')
    fileMsq = open(fnameMsq, 'rb')

    Estore = np.load(fileE)
    Mstore = np.load(fileM)
    Esqstore = np.load(fileEsq)
    Msqstore = np.load(fileMsq)

    # Skiping line with order number & toplogical number  
    next(file)
    
    # Tolopogy number
    c = int(file.readline())
    
    # Going through the cluster files again to read in the 
    # subcluster information
    EOF = False
    while EOF == False:

      line = file.readline().split()

      # Skiping the bond info  
      for b in range(int(line[0])):
        next(file)
        
      # Subgraph info
      sc_Num = int(file.readline())
      for sg in range(sc_Num):
        line = file.readline().split()
        subclusterSize = int(line[0])
        scMult = int(line[1])
        sb_topN = int(line[2])
        
        # In computing contributions from clusters, we first
        # subtract the subcluster weights, except for the single 
        # site subcluster
        for i in range(4):
            weights[i][N][c][:] -= weights[i][subclusterSize][sb_topN][:] * scMult
 
      # We then add the properties computed in ED and subtract the 
      # single site contributions. See https://arxiv.org/pdf/1207.3366.pdf
      # for more details.
      weights[0,N,c,:] += Estore[c,:]   - N*singleE[:]
      weights[1,N,c,:] += Mstore[c,:]   - N*singleM[:]
      weights[2,N,c,:] += Esqstore[c,:] - N*singleEsq[:]
      weights[3,N,c,:] += Msqstore[c,:] - N*singleMsq[:]

      # We are now ready to put together the partial sums, using 
      # the cluster contributions and corresponding lattice constants
      O[:,N,:] += multi[N, c]*weights[:,N, c, :]

      next(file)
      if "Topo. # : Multp." in file.readline():
        # It's the end of the file 
        EOF = True
        file.close()
        fileE.close()
        fileM.close()
        fileEsq.close()
        fileMsq.close()
        break
      else:
        c += 1
  return O

 
O = NLCE()

#---------------------------------------------------------
# In this part, we add partial sums to obtain the NLCE sums
# upto different orders. We can do it without any tricks 
# (raw sums), or using numerical resummation tricks (here
# we have implemented the Wynn algorithm, see Sec. 4 of 
# https://arxiv.org/pdf/1207.3366.pdf for more details).
# The Euler resummation algorithm is also useful and worth
# implementing.

# Doing the raw sums for NLCE
raw_O = np.zeros([4,Nmax+1,NT])
for N in range(1,Nmax+1):
    raw_O[:,N,:] = raw_O[:,N-1,:] + O[:,N,:]

# Doing the sums according to the Wynn algorithm for faster convergence
nwynn = 2
epsln_O = np.zeros([4,Nmax+1,NT,8])
for i in range(NT):
    for N in range(1,Nmax+1):
        epsln_O[:,N,i,0] = epsln_O[:,N-1,i,0] + O[:,N,i]

    for k in range(1,2*nwynn+1):
        for N in range(1,Nmax-k+1):
            for j in range(4):
                delta = epsln_O[j,N+1,i,k-1] - epsln_O[j,N,i,k-1]
                
                #if abs(delta/epsln_O[j,N,i,k]) > 0.0:
                if k == 1:
                    epsln_O[j,N,i,k] = 1.0/delta
                else:
                    epsln_O[j,N,i,k] = epsln_O[j,N+1,i,k-2] + 1.0/delta

nwynn *= 2

plt.figure(figsize=(10,8))
plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('legend', fontsize=15) 

colors = ['g', 'r', 'b', 'c', 'k', 'm', 'y']

# plt.rcParams['font.size']= 15

# This is to optionally skip very low temperatures data in plotting
# as they can get wild. Tlimit sets the index of the lowest 
# temperature ploted.
Tlimit = 20

# Energy
ax = plt.subplot(221)
for N in range(Nmax-nwynn,Nmax-nwynn+1):
    plt.plot(Temp[Tlimit:],epsln_O[0,N,Tlimit:,nwynn],'--',label='Wynn')
                
for N in range(Nmax-1,Nmax+1):
    plt.plot(Temp[Tlimit:],raw_O[0,N,Tlimit:], colors[N-3], label='%s' % N)

plt.ylabel("E")
plt.xlabel("T")
plt.xlim(0.5,100)
plt.ylim(-1.5,0)
plt.xscale('log')
ax.legend(loc='lower right', frameon=False)


# Magnetization
plt.subplot(222)
for N in range(Nmax-nwynn,Nmax-nwynn+1):
    plt.plot(Temp[Tlimit:],epsln_O[1,N,Tlimit:,nwynn],'--')
                
for N in range(Nmax-1,Nmax+1):
    plt.plot(Temp[Tlimit:],raw_O[1,N,Tlimit:], colors[N-3])
plt.xlabel("T")
plt.ylabel("M")
plt.xlim(0.5,100)
plt.ylim(-0.1,0.1)
plt.xscale('log')

    
# Heat capacity
# Compare this quantity with Fig. 14 of https://arxiv.org/pdf/0706.3254.pdf
plt.subplot(223)
for N in range(Nmax-nwynn,Nmax-nwynn+1):
    plt.plot(Temp[Tlimit:],(epsln_O[2,N,Tlimit:,nwynn]-epsln_O[0,N,Tlimit:,nwynn]**2)/Temp[Tlimit:]**2,'--')
                
for N in range(Nmax-1,Nmax+1):
    plt.plot(Temp[Tlimit:],raw_O[2,N,Tlimit:]/Temp[Tlimit:]**2, colors[N-3])
plt.xlabel("T")
plt.ylabel(r'$C_{v}$')
plt.xlim(0.5,100)
plt.ylim(0,0.2)
plt.xscale('log')
    
    
# Magnetic susceptibility
plt.subplot(224)
for N in range(Nmax-nwynn,Nmax-nwynn+1):
    plt.plot(Temp[Tlimit:],epsln_O[3,N,Tlimit:,nwynn]/Temp[Tlimit:],'--')
                
for N in range(Nmax-1,Nmax+1):
    plt.plot(Temp[Tlimit:],raw_O[3,N,Tlimit:]/Temp[Tlimit:], colors[N-3])
plt.xlim(0.5,100)
plt.ylim(-1,1)
plt.xlabel("T")
plt.ylabel(r'$X_{m}$', labelpad=-11)
plt.xscale('log')
#plt.yscale('log')

plt.legend(loc='lower right')
plt.tight_layout()
#plt.subplots_adjust(right = 1.1)
plt.show()





