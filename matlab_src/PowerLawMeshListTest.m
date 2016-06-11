ccc

matlab_config_filename   = '../config/ConfigurationMatlab.cfg';
cfg_template_list = 'june10_templatelist';
Nmeshes   = 24; % number of meshes to be generated

in = fopen(cfg_template_list,'r');
config_template_filename = fgetl(in);

while (config_template_filename~=-1)
    
    [path,runname,ext] = fileparts(config_template_filename);
    config_list_filename = [runname '_runlist'];
    
    Files.matlab_config_filename   = matlab_config_filename;
    Files.config_template_filename = config_template_filename;
    Files.config_list_filename     = config_list_filename;
    
    % Ceres topography spectrum power law parameters
    
    cfg = ReadConfig(Files);
    
    tic
    GeneratePowerLawMesh2(Files,cfg,Nmeshes,runname);
    toc
   
    config_template_filename = fgetl(in);
    
    close all
end

