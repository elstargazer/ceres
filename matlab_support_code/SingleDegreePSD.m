function lmcosi_shape = SingleDegreePSD(r_mean,beta,intercept,L,n_only)

% r_mean = 470000;
% beta = -3.5;
% intercept = 8.9;
% L = 100;
% n_only = 6;

lmcosi_shape = CreateEmptylmcosi(L);
lmcosi_shape(1,3) = r_mean;
lmcosi_shape(2,3:4) = 0;
lmcosi_shape(3,3:4) = 0;

for n=2:2:lmcosi_shape(end,1)
    rand_coeff = ((rand<0.5)*2-1);
    %force the deg 2 and 4 coefficients to be -1 and +1 to make body
    %rotationally stable and avoid unrealistic mountain at pole
    if(n==2)
        rand_coeff = -1;
    end
    if(n==4)
        rand_coeff = 1;
    end
    
    %         fileID = fopen(info_file_name, 'a');
    %         formatSpec = '%2.1f ';
    %         fprintf(fileID, formatSpec, rand_coeff);
    %         ii=1:1;
    %         fprintf(fileID, '\n', ii);
    %         fclose(fileID);
    
    if n==n_only
        only_deg = 1;
    else
        only_deg = 0;
    end
    
    lmcosi_shape((n+1)*n/2+1,3) = only_deg*rand_coeff*sqrt((2*n+1)*...
        10^polyval([beta intercept],log10(n)));
end