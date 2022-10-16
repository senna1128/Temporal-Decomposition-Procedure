# This function computes the KKT residual

function ComputeKKT(N,C1,C2,dk,xuout,lambdaout)
    GradXU, GradF = zeros(2*N+1), zeros(N+1)
    for i = 1:N
        GradXU[2*i-1]=2*C1*(xuout[2*i-1]-dk[i])-2*sin(2*(xuout[2*i-1]-dk[i]))+lambdaout[i]-lambdaout[i+1]
        GradXU[2*i] = -2*C2*(xuout[2*i]-dk[i]) - lambdaout[i+1]
        GradF[i+1] = xuout[2*i+1]-(xuout[2*i-1]+xuout[2*i]+dk[i])
    end
    GradF[1] = xuout[1]
    GradXU[end] = 2*C1*xuout[end]+lambdaout[end]
    return GradXU, GradF, norm([GradXU;GradF])
end
