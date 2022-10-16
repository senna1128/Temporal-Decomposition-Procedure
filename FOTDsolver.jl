## This fucntion conducts one step Newton

function FOTDsolver(N,B,T,TempDecom,OverTempDecom,DiagH,subG,GradXU,GradF,mu,id_Warm=1)
    ## Define Result Vector
    XU_dir, Lam_dir = zeros(2*N+1), zeros(N+1)

    if id_Warm == 1
        ## First subproblem
        m1, m2 = 1, 2*OverTempDecom[1,2]-1
        n1, n2 = 1, TempDecom[1,2]-1
        H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
        H_sub[end]+= mu
        GradF_sub = GradF[1:OverTempDecom[1,2]]
        subG_sub = subG[1:OverTempDecom[1,2],1:m2]
        HH = Diagonal(H_sub)
        KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(OverTempDecom[1,2],OverTempDecom[1,2])))
        Grad_sub = vcat(GradXU_sub,GradF_sub)
        Itera = lu(KK)\-Grad_sub
        XU_dir[1:2*n2] = Itera[1:2*n2]
        Lam_dir[1:n2] = Itera[m2+1:m2+n2]
        ## Middle subproblem
        for i=2:T-1
            m1, m2 = 2*OverTempDecom[i,1]-1, 2*OverTempDecom[i,2]-1
            n1, n2 = TempDecom[i,1], TempDecom[i,2]-1
            sub_len = OverTempDecom[i,2]-OverTempDecom[i,1]+1
            H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
            H_sub[end]+= mu
            subG_sub = subG[1:sub_len,1:2*sub_len-1]
            GradF_sub = GradF[OverTempDecom[i,1]:OverTempDecom[i,2]]
            HH = Diagonal(H_sub)
            KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(sub_len,sub_len)))
            Grad_sub = vcat(GradXU_sub,GradF_sub)
            Itera = lu(KK)\-Grad_sub
            XU_dir[2*n1-1:2*n2] = Itera[2*B+1:2*(n2-n1+B+1)]
            Lam_dir[n1:n2] = Itera[2*sub_len+B:end-B-1]
        end
        ## Last subproblem
        m1, m2 = 2*OverTempDecom[end,1]-1, 2*N+1
        n1, n2 = TempDecom[end,1], N+1
        sub_len = N-OverTempDecom[end,1]+2
        H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
        subG_sub = subG[1:sub_len,1:2*sub_len-1]
        GradF_sub = GradF[OverTempDecom[end,1]:n2]
        HH = Diagonal(H_sub)
        KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(sub_len,sub_len)))
        Grad_sub = vcat(GradXU_sub,GradF_sub)
        Itera = lu(KK)\-Grad_sub
        XU_dir[2*n1-1:2*n2-1] = Itera[2*B+1:2*(B+n2-n1)+1]
        Lam_dir[n1:n2] = Itera[2*sub_len+B:end]
    elseif id_Warm == 2
        ## First subproblem
        m1, m2 = 1, 2*OverTempDecom[1,2]-1
        n1, n2 = 1, TempDecom[1,2]-1
        H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
        H_sub[end]+= mu
        GradF_sub = GradF[1:OverTempDecom[1,2]]
        subG_sub = subG[1:OverTempDecom[1,2],1:m2]
        HH = Diagonal(H_sub)
        KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(OverTempDecom[1,2],OverTempDecom[1,2])))
        Grad_sub = vcat(GradXU_sub,GradF_sub)
        Itera = gmres(KK,-Grad_sub)
        XU_dir[1:2*n2] = Itera[1:2*n2]
        Lam_dir[1:n2] = Itera[m2+1:m2+n2]
        ## Middle subproblem
        for i=2:T-1
            m1, m2 = 2*OverTempDecom[i,1]-1, 2*OverTempDecom[i,2]-1
            n1, n2 = TempDecom[i,1], TempDecom[i,2]-1
            sub_len = OverTempDecom[i,2]-OverTempDecom[i,1]+1
            H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
            H_sub[end]+= mu
            subG_sub = subG[1:sub_len,1:2*sub_len-1]
            GradF_sub = GradF[OverTempDecom[i,1]:OverTempDecom[i,2]]
            HH = Diagonal(H_sub)
            KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(sub_len,sub_len)))
            Grad_sub = vcat(GradXU_sub,GradF_sub)
            Itera = gmres(KK,-Grad_sub)
            XU_dir[2*n1-1:2*n2] = Itera[2*B+1:2*(n2-n1+B+1)]
            Lam_dir[n1:n2] = Itera[2*sub_len+B:end-B-1]
        end
        ## Last subproblem
        m1, m2 = 2*OverTempDecom[end,1]-1, 2*N+1
        n1, n2 = TempDecom[end,1], N+1
        sub_len = N-OverTempDecom[end,1]+2
        H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
        subG_sub = subG[1:sub_len,1:2*sub_len-1]
        GradF_sub = GradF[OverTempDecom[end,1]:n2]
        HH = Diagonal(H_sub)
        KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(sub_len,sub_len)))
        Grad_sub = vcat(GradXU_sub,GradF_sub)
        Itera = gmres(KK,-Grad_sub)
        XU_dir[2*n1-1:2*n2-1] = Itera[2*B+1:2*(B+n2-n1)+1]
        Lam_dir[n1:n2] = Itera[2*sub_len+B:end]
    else
        ## First subproblem
        m1, m2 = 1, 2*OverTempDecom[1,2]-1
        n1, n2 = 1, TempDecom[1,2]-1
        H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
        H_sub[end]+= mu
        GradF_sub = GradF[1:OverTempDecom[1,2]]
        subG_sub = subG[1:OverTempDecom[1,2],1:m2]
        HH = Diagonal(H_sub)
        KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(OverTempDecom[1,2],OverTempDecom[1,2])))
        Grad_sub = vcat(GradXU_sub,GradF_sub)
        Itera = idrs(KK,-Grad_sub)
        XU_dir[1:2*n2] = Itera[1:2*n2]
        Lam_dir[1:n2] = Itera[m2+1:m2+n2]
        ## Middle subproblem
        for i=2:T-1
            m1, m2 = 2*OverTempDecom[i,1]-1, 2*OverTempDecom[i,2]-1
            n1, n2 = TempDecom[i,1], TempDecom[i,2]-1
            sub_len = OverTempDecom[i,2]-OverTempDecom[i,1]+1
            H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
            H_sub[end]+= mu
            subG_sub = subG[1:sub_len,1:2*sub_len-1]
            GradF_sub = GradF[OverTempDecom[i,1]:OverTempDecom[i,2]]
            HH = Diagonal(H_sub)
            KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(sub_len,sub_len)))
            Grad_sub = vcat(GradXU_sub,GradF_sub)
            Itera = idrs(KK,-Grad_sub)
            XU_dir[2*n1-1:2*n2] = Itera[2*B+1:2*(n2-n1+B+1)]
            Lam_dir[n1:n2] = Itera[2*sub_len+B:end-B-1]
        end
        ## Last subproblem
        m1, m2 = 2*OverTempDecom[end,1]-1, 2*N+1
        n1, n2 = TempDecom[end,1], N+1
        sub_len = N-OverTempDecom[end,1]+2
        H_sub, GradXU_sub = DiagH[m1:m2], GradXU[m1:m2]
        subG_sub = subG[1:sub_len,1:2*sub_len-1]
        GradF_sub = GradF[OverTempDecom[end,1]:n2]
        HH = Diagonal(H_sub)
        KK = hcat(vcat(HH,subG_sub),vcat(transpose(subG_sub),zeros(sub_len,sub_len)))
        Grad_sub = vcat(GradXU_sub,GradF_sub)
        Itera = idrs(KK,-Grad_sub)
        XU_dir[2*n1-1:2*n2-1] = Itera[2*B+1:2*(B+n2-n1)+1]
        Lam_dir[n1:n2] = Itera[2*sub_len+B:end]
    end
    return XU_dir, Lam_dir
end
