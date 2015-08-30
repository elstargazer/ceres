ccc
matlab_config_filename   = '~/Dawn/FE/config/ConfigurationMatlab.cfg';
config_template_filename = '~/Dawn/FE/config/ConfigurationTemplate.cfg';
Nmeshes   = 10; % number of meshes to be generated

% Ceres topography spectrum power law parameters
r_mean    = 470000;
beta      = -3.72;
intercept = 8.079;
GeneratePowerLawMesh(matlab_config_filename,config_template_filename,...
    r_mean,beta,intercept,Nmeshes);