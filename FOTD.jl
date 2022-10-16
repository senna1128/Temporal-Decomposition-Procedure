include("FOTDNLP.jl")
## This fucntion use FOTD to solve nonlinear program

function FOTD(N,C1,C2,dk,MM,BB,FOTD_mu,FOTD_eta1,FOTD_eta2,FOTD_beta,FOTD_rho,Rep,XU_Init,Lam_Init)
    LenMM, LenBB, LenMu = length(MM), length(BB), length(FOTD_mu)
    ## Result Vectors
    FOTDKKT = reshape([[] for i = 1:LenMM*LenBB*LenMu*3],(LenMM,LenBB,LenMu,3))
    FOTDTime = reshape([[] for i = 1:LenMM*LenBB*LenMu*3],(LenMM,LenBB,LenMu,3))
    FOTDnotconv = reshape([0 for i = 1:LenMM*LenBB*LenMu*3],(LenMM,LenBB,LenMu,3))
    ## Go over all cases
    for id_MM = 1:LenMM
        for id_BB = 1:LenBB
            for id_Mu = 1:LenMu
                for id_Warm = 1:3
                    for rrep = 1:Rep
                        FOTDKKTVec,FOTDTimeVec,IdConv = FOTDNLP(N,C1,C2,dk,MM[id_MM],BB[id_BB],FOTD_mu[id_Mu],FOTD_eta1,FOTD_eta2,FOTD_beta,FOTD_rho,XU_Init[:,rrep],Lam_Init[:,rrep],id_Warm)
                        if IdConv == 1
                            push!(FOTDKKT[id_MM,id_BB,id_Mu,id_Warm],FOTDKKTVec)
                            push!(FOTDTime[id_MM,id_BB,id_Mu,id_Warm],FOTDTimeVec)
                        else
                            FOTDnotconv[id_MM,id_BB,id_Mu,id_Warm] += 1
                        end
                    end
                end
            end
        end
    end
    # return results
    return FOTDKKT, FOTDTime, FOTDnotconv
end
