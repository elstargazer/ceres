function config = ReadConfig(Files)

in = fopen(Files.config_template_filename,'r');

str = 1;

key_r_mean = 'r_mean';
key_period = 'period';
key_beta   = 'beta';
key_intercept = 'intercept';

key_depths_eta = 'depths_eta';
key_eta_kinks = 'eta_kinks';
key_depths_rho = 'depths_rho';
key_rho = 'rho';
key_mat_id = 'material_id';


max_elems = 10;

while (str~=-1)
    str = strtrim(fgetl(in));
    
    if strncmpi(str,key_period,numel(key_period))
        
        config.T = sscanf(str,'%*s = %f;');
        
    elseif strncmpi(str,key_r_mean,numel(key_r_mean))
        
        config.r_mean = sscanf(str,'%*s = %f;');
        
    elseif strncmpi(str,key_beta,numel(key_beta))
        
        config.beta = sscanf(str,'%*s = %f;');
        
    elseif strncmpi(str,key_intercept,numel(key_intercept))
        
        config.intercept = sscanf(str,'%*s = %f;');
        
    elseif strncmpi(str,key_depths_eta,numel(key_depths_eta))
        
        config.depths_eta = sscanf(str,['%*s = [' repmat('%f, ',[1 max_elems]) ' %f];']);
        
    elseif strncmpi(str,key_eta_kinks,numel(key_eta_kinks))
        
        config.eta_kinks = sscanf(str,['%*s = [' repmat('%E,',[1 max_elems]) '%E];']);
    
    elseif strncmpi(str,key_depths_rho,numel(key_depths_rho))
        
        config.depths_rho = sscanf(str,['%*s = [' repmat('%f, ',[1 max_elems]) ' %f];']);
        
    elseif strncmpi(str,key_rho,numel(key_rho))
        
        config.rho = sscanf(str,['%*s = [' repmat('%f, ',[1 max_elems]) '%f];']);
        
    elseif strncmpi(str,key_mat_id,numel(key_mat_id))
        
        config.mat_id = sscanf(str,['%*s = [' repmat('%d, ',[1 max_elems]) '%f];']);
        
        
    else
        
        
        
    end
end



