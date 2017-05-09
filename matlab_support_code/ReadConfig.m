function config = ReadConfig(Files)

in = fopen(Files.config_template_filename,'r');

str = 1;

key_r_mean = 'r_mean';
key_period = 'period';
key_beta   = 'beta';
key_intercept = 'intercept';
spectrum_filename = 'spectrum_filename';

key_depths_eta = 'depths_eta';
key_eta_kinks = 'eta_kinks';
key_depths_rho = 'depths_rho';
key_rho = 'rho';
key_mat_id = 'material_id';
key_cell_height = 'cell_height';
key_output = 'output_folder';


max_elems = 10;

str = fgetl(in);

while (str~=-1)
    
    str = strtrim(str);
    
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
        
    elseif strncmpi(str,key_cell_height,numel(key_cell_height))
        
        config.cell_height = sscanf(str,['%*s = [' repmat('%f, ',[1 max_elems]) ' %f];']);
        
    elseif strncmpi(str,key_eta_kinks,numel(key_eta_kinks))
        
        config.eta_kinks = sscanf(str,['%*s = [' repmat('%E,',[1 max_elems]) '%E];']);
        
    elseif strncmpi(str,key_depths_rho,numel(key_depths_rho))
        
        config.depths_rho = sscanf(str,['%*s = [' repmat('%f, ',[1 max_elems]) ' %f];']);
        
    elseif strncmpi(str,key_rho,numel(key_rho))
        
        config.rho = sscanf(str,['%*s = [' repmat('%f, ',[1 max_elems]) '%f];']);
        
    elseif strncmpi(str,key_mat_id,numel(key_mat_id))
        
        config.mat_id = sscanf(str,['%*s = [' repmat('%d, ',[1 max_elems]) '%f];']);
        
    elseif strncmpi(str,spectrum_filename,numel(spectrum_filename))
        
        config.spectrum_filename = sscanf(str,'%*s = "%s"');
        config.spectrum_filename = config.spectrum_filename(1:end-2);
        
    elseif strncmpi(str,key_output,numel(key_mat_id))
        
        textline = sscanf(str,'%*s = %s;');
        subchunk = regexp(textline, '(?<=")[^"]+(?=")', 'match');
        config.output_folder = subchunk{1};
        
        
    else
        
        
        
    end
    
    str = fgetl(in);
    
end

fclose(in);



