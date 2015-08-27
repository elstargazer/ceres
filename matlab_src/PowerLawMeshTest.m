% ccc
config_filename = '~/Dawn/FE/config/ConfigurationMatlab.cfg';

Nmeshes   = 1; % number of meshes to be generated

% Ceres topography spectrum power law parameters
r_mean    = 470000;
beta      = -3.72;
intercept = 8.079;
GeneratePowerLawMesh(config_filename,r_mean,beta,intercept,2)