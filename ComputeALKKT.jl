# This function computes the KKT and Gradient of Augmented Lagrangian

function ComputeALKKT(N,C1,C2,dk,FOTD_eta1,FOTD_eta2,XU,Lam)
    GradXU,GradF,DiagH = zeros(2*N+1),zeros(N+1),zeros(2*N+1)
    HnabZ,GTf,GnabZ = zeros(2*N+1),zeros(2*N+1),zeros(N+1)
    ## Initial Stage
    GradXU[1] = 2*C1*(XU[1]-dk[1])-2*sin(2*(XU[1]-dk[1]))+Lam[1]-Lam[2]
    GradXU[2] = -2*C2*(XU[2]-dk[1])-Lam[2]
    GradF[1],GradF[2] = XU[1], XU[3]-(XU[1]+XU[2]+dk[1])
    HnabZ[1],HnabZ[2] = (2*C1-4*cos(2*(XU[1]-dk[1])))*GradXU[1],-2*C2*GradXU[2]
    DiagH[1],DiagH[2] = 2*C1-4*cos(2*(XU[1]-dk[1])),-2*C2
    GTf[1],GTf[2] = GradF[1]-GradF[2],-GradF[2]
    GnabZ[1] = GradXU[1]
    ## Middle Stage
    for i = 2:N
        GradXU[2*i-1]=2*C1*(XU[2*i-1]-dk[i])-2*sin(2*(XU[2*i-1]-dk[i]))+Lam[i]-Lam[i+1]
        GradXU[2*i] = -2*C2*(XU[2*i]-dk[i])-Lam[i+1]
        GradF[i+1] = XU[2*i+1]-(XU[2*i-1]+XU[2*i]+dk[i])
        HnabZ[2*i-1] = (2*C1-4*cos(2*(XU[2*i-1]-dk[i])))*GradXU[2*i-1]
        HnabZ[2*i] = -2*C2*GradXU[2*i]
        DiagH[2*i-1],DiagH[2*i] = 2*C1-4*cos(2*(XU[2*i-1]-dk[i])),-2*C2
        GTf[2*i-1],GTf[2*i] = GradF[i]-GradF[i+1],-GradF[i+1]
        GnabZ[i] = GradXU[2*i-1]-(GradXU[2*i-2]+GradXU[2*i-3])
    end
    GradXU[end] = 2*C1*XU[end]+Lam[end]
    HnabZ[end],DiagH[end],GTf[end] = 2*C1*GradXU[end],2*C1,GradF[end]
    GnabZ[end] = GradXU[end]-(GradXU[end-1]+GradXU[end-2])

    return GradXU,GradF,GradXU+FOTD_eta2*HnabZ+FOTD_eta1*GTf,GradF+FOTD_eta2*GnabZ,DiagH
end
