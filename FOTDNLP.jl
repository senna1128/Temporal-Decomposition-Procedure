include("ConsSubG.jl")
include("FOTDsolver.jl")
include("backtrack.jl")
include("ComputeALKKT.jl")

function FOTDNLP(N,C1,C2,dk,M,B,mu,FOTD_eta1,FOTD_eta2,FOTD_beta,FOTD_rho,XU,Lam,id_Warm=1,IterMax=40,epsthres=1e-6)
    # compute number of subproblems
    T = Int(floor(N/M))
    # temporal decomposition
    TempDecom,OverTempDecom = Array{Int}(undef,T,2),Array{Int}(undef,T,2)
    for i = 1:T
        TempDecom[i,:],OverTempDecom[i,:] = [(i-1)*M,i*M],[max((i-1)*M-B,1), min(i*M+B,N)]
    end
    TempDecom[1,1],TempDecom[end,end] = 1,N+1
    OverTempDecom[1,1],OverTempDecom[end,end] = 1,N+1
    # Result Vector
    FOTDKKTVec,IdConv,subG = [],0,ConsSubG(2*(M+B)+1)

    ## Compute Augmented Lagrangian Gradient etc.
    GradXU,GradF,GradALXU,GradALF,DiagH = ComputeALKKT(N,C1,C2,dk,FOTD_eta1,FOTD_eta2,XU,Lam)
    push!(FOTDKKTVec,norm([GradXU;GradF]))
    eps,tau = FOTDKKTVec[end], 1
    # start the iteration
    start = time()
    while tau<=IterMax && eps>epsthres
        # Obtain search direction
        XU_dir,Lam_dir = FOTDsolver(N,B,T,TempDecom,OverTempDecom,DiagH,subG,GradXU,GradF,mu,id_Warm)
        # Backtracking line search
        step = backtrack(N,C1,C2,dk,XU,Lam,XU_dir,Lam_dir,GradXU,GradF,GradALXU,GradALF,FOTD_eta1,FOTD_eta2,FOTD_beta,FOTD_rho)
        XU = XU + step*XU_dir
        Lam = Lam + step*Lam_dir
        ## Compute Augmented Lagrangian Gradient etc.
        GradXU,GradF,GradALXU,GradALF,DiagH = ComputeALKKT(N,C1,C2,dk,FOTD_eta1,FOTD_eta2,XU,Lam)
        push!(FOTDKKTVec,norm([GradXU;GradF]))
        eps,tau = step*norm([XU_dir;Lam_dir]),tau+1
        if eps <= epsthres || FOTDKKTVec[end]<=epsthres
            IdConv = 1
            break
        end
    end
    opttime = time()-start
    if IdConv == 1
        return FOTDKKTVec, opttime, 1
    else
        return [],0,0
    end
end
