
function ConsSubG(MMBB)
    subG = zeros(MMBB,2*MMBB-1)
    subG[1,1]=1
    for i = 2:MMBB
        subG[i,2*i-1],subG[i,2*i-2],subG[i,2*i-3] = 1,-1,-1
    end
    return subG
end
