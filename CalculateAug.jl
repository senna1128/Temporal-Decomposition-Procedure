## This fucntion calculates the augmented Lagrange function value

function CalculateAug(N,C1,C2,dk,XU,Lam,GradXU,GradF,FOTD_eta1,FOTD_eta2,Id_new)
    if Id_new == 0
        Obj1 = sum(2*cos(XU[2*i-1]-dk[i])^2+C1*(XU[2*i-1]-dk[i])^2-C2*(XU[2*i]-dk[i])^2 for i = 1:N)+C1*XU[end]^2
        Obj2 = Lam'GradF
        Obj3 = FOTD_eta1/2*norm(GradF)^2
        Obj4 = FOTD_eta2/2*norm(GradXU)^2
        return Obj1+Obj2+Obj3+Obj4
    else
        Obj1 = sum(2*cos(XU[2*i-1]-dk[i])^2+C1*(XU[2*i-1]-dk[i])^2-C2*(XU[2*i]-dk[i])^2 for i = 1:N)+C1*XU[end]^2
        gradXU, gradF,~ = ComputeKKT(N,C1,C2,dk,XU,Lam)
        Obj2 = Lam'gradF
        Obj3 = FOTD_eta1/2*norm(gradF)^2
        Obj4 = FOTD_eta2/2*norm(gradXU)^2
        return Obj1+Obj2+Obj3+Obj4
    end
end
