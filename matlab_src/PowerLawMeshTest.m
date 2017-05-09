ccc

matlab_config_filename   = '/Users/ermakov/Dawn/FE/config/ConfigurationMatlab.cfg';
config_template_filename = '/Users/ermakov/Dawn/FE/config/May5_g10_softint.cfg';
config_list_filename     = '/Users/ermakov/Dawn/FE/May5_g10_softint_runlist';
runname                  = 'May5_g10_softint';

Files.matlab_config_filename   = matlab_config_filename;
Files.config_template_filename = config_template_filename;
Files.config_list_filename     = config_list_filename;

Nmeshes   = 24; % number of meshes to be generated

% Ceres topography spectrum power law parameters

cfg = ReadConfig(Files);

tic
GeneratePowerLawMesh2(Files,cfg,Nmeshes,runname);
toc
