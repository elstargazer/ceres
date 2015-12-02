function lmcosi_shape = PowerLawSH(r_mean,beta,intercept,L)

lmcosi_shape = CreateEmptylmcosi(L);
lmcosi_shape(1,3) = r_mean;
lmcosi_shape(2,3:4) = 0;
lmcosi_shape(3,3:4) = 0;

for n=2:2:lmcosi_shape(end,1)
    lmcosi_shape((n+1)*n/2+1,3) = ((rand<0.5)*2-1)*sqrt((2*n+1)*...
        10^polyval([beta intercept],log10(n)));
end