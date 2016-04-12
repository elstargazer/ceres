ccc

matlab_config_filename   = '../config/ConfigurationMatlab.cfg';
config_template_filename = '../config/Run115.cfg';
config_list_filename     = '../RunList.txt';
runname                  = 'run115';

Files.matlab_config_filename   = matlab_config_filename;
Files.config_template_filename = config_template_filename;
Files.config_list_filename     = config_list_filename;

Nmeshes   = 24; % number of meshes to be generated

% Ceres topography spectrum power law parameters

cfg = ReadConfig(Files);

GeneratePowerLawMesh2(Files,cfg,Nmeshes,runname);
