include("CalculateAug.jl")
## This fucntion implements the backtracking line search method

function backtrack(N,C1,C2,dk,XU,Lam,XU_dir,Lam_dir,GradXU,GradF,GradALXU,GradALF,FOTD_eta1,FOTD_eta2,FOTD_beta,FOTD_rho,thres=1e-5)
    quant1 = CalculateAug(N,C1,C2,dk,XU,Lam,GradXU,GradF,FOTD_eta1,FOTD_eta2,0)
    quant2 = FOTD_beta*(GradALXU'XU_dir+GradALF'Lam_dir)
    step = 1
    while CalculateAug(N,C1,C2,dk,XU+step*XU_dir,Lam+step*Lam_dir,0,0,FOTD_eta1,FOTD_eta2,1)>quant1+step*quant2 && step>thres
        step *= FOTD_rho
    end
    return step
end
