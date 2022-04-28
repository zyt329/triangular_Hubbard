import julia
from julia import Main
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from itertools import permutations

# In this part, we read in the bond information for every
# cluster and solve the Ising model for them.

# Maximum NLCE order
Nmax = 8

# hopping t
Main.t = 1.0

# Interaction U
Main.U = 8.0

# Holding "lattice constants" (multiplicities) for all clusters
# in all orders. Probably good for up to order 9. Need to check
# what is needed exactly for the size of the second dimension.
multi = np.zeros([Nmax + 1, 3500], int)

# Holding the two site numbers and optionally the bond type
# for every bond.
site1 = np.zeros(40, int)
site2 = np.zeros(40, int)
typ = np.zeros(40, int)

# Exchange constant and the magnetic field
# J = np.array([1, 0])
# Bfield = 0.0


# Counter for isomorphically distinct clusters
topN = 0

# load julia script
Main.include("./Hubbard.jl")
Main.include("./measurements_Hubbard.jl")

Main.name = "order=" + str(Nmax) + "_U=" + str(Main.U)

# Loop over the NLCE order
for N in range(2, Nmax + 1):
    # Passing N to the julia script
    Main.N = N

    # Generate sector info for julia script
    Main.eval("sectors_info = sectors_info_gen(N=N)")

    # Change 1 to 2 in the line below and in the next cell
    # to include n.n.n. bonds
    floc = "./NLCE_Clusters/"
    fbase = floc + "NLCE_1_" + str(N)
    fname = fbase + ".txt"
    file = open(fname, "r")

    # Skips line with order number
    next(file)

    # Going through the file for each order line by line
    # and reading in the cluster information
    topN = int(file.readline())
    print("ORDER", N)

    EOF = False
    while EOF == False:
        # passing topN to julia script
        Main.NTOP = topN

        line = file.readline().split()

        # Get the number of bonds from this line, and pass it to julia script
        nB = int(line[0])
        Main.nB = nB

        # Read in the bond information, and pass it to julia script
        Main.eval("site1 = Int64[]")
        Main.eval("site2 = Int64[]")
        for b in range(nB):
            line = file.readline().split()
            site1[b] = int(eval(line[0]))
            Main.eval(f"push!(site1, {site1[b]} + 1)")
            site2[b] = int(eval(line[1]))
            Main.eval(f"push!(site2, {site2[b]} + 1)")
            typ[b] = int(eval(line[2])) - 1

        Main.eval("bonds = [[site1[i], site2[i]] for i in 1:nB]")
        # calculate eigenvalues of matrix and expectations of quantities in eigenstates
        Main.eval(
            "quantities = E_Quants(N=N, t=t, U=U, sectors_info=sectors_info, bonds=bonds)"
        )
        # printing calculated quantities to a file for save
        Main.eval(
            "printing(quantities; Quant_names = quant_names, name = name, NTOP=NTOP, N=N)"
        )

        # ---------------------------------------------------------

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
            # next(file)
            for i in range(topN + 1):
                multi[N][i] = file.readline().split()[1]
            EOF = True
            file.close()
            break
        else:
            # If not yet done with clusters
            topN += 1

