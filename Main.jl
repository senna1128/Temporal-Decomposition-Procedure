
#### Main File For Simulation
## Load packages
using NLPModels
using JuMP
using LinearOperators
using OptimizationProblems
using MathProgBase
using ForwardDiff
using NLPModelsJuMP
using LinearAlgebra
using Distributed
using Ipopt
using DataFrames
using PyPlot
using MATLAB
using ADNLPModels
using Glob
using DelimitedFiles
using Random
using Distributions
using NLPModelsIpopt
using Statistics
using IterativeSolvers

## Change working directory
cd("$(homedir())/...")

# Pamareter class
module Parameter
    struct ProbParams
        verbose                # Do we create dump dir?
        solver                 # Ipopt, PipsNlp
        rep::Int               # replication times
        N::Int                 # horizon length
        MM::Array{Int}         # short horizon length
        BB::Array{Int}         # overlap size
        C1::Float64            # coefficient 1
        C2::Float64            # coefficient 2, 4C2<C1-2
        dk::Array              # reference variable dk, N by 1 vector
    end
    struct FOTDParams
        mu::Array{Float64}      # scale Parameter
        eta1::Float64           # augmented Lagrange penalty parameter
        eta2::Float64           # augmented Lagrange penalty parameter
        beta::Float64           # armijo parameter between 0 and 0.5
        rho::Float64            # backtracking line search parameter between 0 and 1
    end
end

using Main.Parameter


include("FOTD.jl")


function _main(ProbParams,FOTDParams)
    ## Problem parameter
    Prob_solver,Rep = ProbParams.solver,ProbParams.rep
    N,MM,BB = ProbParams.N,ProbParams.MM,ProbParams.BB
    C1,C2,dk = ProbParams.C1, ProbParams.C2, ProbParams.dk
    ## Initialization
    XU_Init,Lam_Init = rand(Uniform(-1e5,1e5),(2*N+1,Rep)),rand(Uniform(-1e5,1e5),(N+1,Rep))
    XU_Init[:,1],Lam_Init[:,1] = zeros(2*N+1),zeros(N+1)
    XU_Init[1,:] = zeros(Rep)
    ## Algorithmic parameter
    FOTD_mu,FOTD_eta1,FOTD_eta2 = FOTDParams.mu,FOTDParams.eta1,FOTDParams.eta2
    FOTD_beta,FOTD_rho =  FOTDParams.beta,FOTDParams.rho

    ## Solve FOTD
    FOTDKKT, FOTDTime, FOTDnotconv = FOTD(N,C1,C2,dk,MM,BB,FOTD_mu,FOTD_eta1,FOTD_eta2,FOTD_beta,FOTD_rho,Rep,XU_Init,Lam_Init)
end

# main function
function main()
    AllParam = glob("....")
    include(AllParam[CaseId])
    _main(ProbParams,FOTDParams)
end

main()
