ccc

matlab_config_filename   = '../config/ConfigurationMatlab.cfg';
config_template_filename = '../config/Vary_core_g10_55km.cfg';
config_list_filename     = '../Vary_core_g10_55km_runlist';
runname                  = 'Vary_core_g10_55km';

Files.matlab_config_filename   = matlab_config_filename;
Files.config_template_filename = config_template_filename;
Files.config_list_filename     = config_list_filename;

Nmeshes   = 24; % number of meshes to be generated

% Ceres topography spectrum power law parameters

cfg = ReadConfig(Files);

tic
GeneratePowerLawMesh2(Files,cfg,Nmeshes,runname);
toc
