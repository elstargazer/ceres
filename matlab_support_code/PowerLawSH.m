function lmcosi_shape = PowerLawSH(varargin)

if (nargin==5)
    
    %% if power law parameters are given
    r_mean = varargin{1};
    beta = varargin{2};
    intercept = varargin{3};
    L = varargin{4};
    info_file_name = varargin{5};
    
    lmcosi_shape = CreateEmptylmcosi(L);
    lmcosi_shape(1,3) = r_mean;
    lmcosi_shape(2,3:4) = 0;
    lmcosi_shape(3,3:4) = 0;
    
%     fileID = fopen(info_file_name, 'w');
%     fclose(fileID);
    
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
        
        lmcosi_shape((n+1)*n/2+1,3) = rand_coeff*sqrt((2*n+1)*...
            10^polyval([beta intercept],log10(n)));
    end
    
    
elseif (nargin==4)
    
    %% if spherctrum itself is given
    
    r_mean = varargin{1};
    spectrum_filename = varargin{2};
    L = varargin{3};
    info_file_name = varargin{4};
    
    psd = load(spectrum_filename);
    n_psd = 0:numel(psd)-1;
    
    lmcosi_shape = CreateEmptylmcosi(L);
    lmcosi_shape(1,3) = r_mean;
    lmcosi_shape(2,3:4) = 0;
    lmcosi_shape(3,3:4) = 0;
    
%     fileID = fopen(info_file_name, 'w');
%     fclose(fileID);
    
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
        
        %         psd_temp = 10^polyval([beta intercept],log10(n));
        psd_temp = psd(n_psd == n);
        
        lmcosi_shape((n+1)*n/2+1,3) = rand_coeff*sqrt((2*n+1)*...
            psd_temp);
    end
 
end