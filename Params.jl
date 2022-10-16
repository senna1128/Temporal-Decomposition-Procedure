

ProbParams = Parameter.ProbParams(true,
                        "Ipopt",            # solver
                        10,                  # rep
                        5000,               # N
                        [50],               # MM
                        [1,5,25],          # BB
                        8,                  # C1
                        1,                  # C2
                        ones(5000))         # dk

FOTDParams = Parameter.FOTDParams(
                        [1,25,125],         # mu
                        10,                 # eta1
                        0.1,                # eta2
                        0.1,                # beta
                        0.9)                # rho
