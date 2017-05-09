ccc

%% Image settings
fntsize = 20;
fntsize_sm = 10;
im_size=[0 0 20 20];
fig_folder='~/Dawn/Figures/';

%% Input parameters

L = 16;
Rref = 470000;
n_only = 8;

%% Files
matlab_config_filename   = '~/Dawn/FE/config/ConfigurationMatlab.cfg';
config_template_filename = '~/Dawn/FE/config/deg6_42km/deg6_42km_1.cfg';
runname                  = 'deg6_42km';
filename_real_spc        = '~/Dawn/FORTRAN/PSD_SPG_HAMO_20160107.txt';

%% Read configuration

Files.matlab_config_filename   = matlab_config_filename;
Files.config_template_filename = config_template_filename;

% read matlab config
in = fopen(matlab_config_filename);
str = fscanf(in,'FE_folder = %s\n',1);
FE_folder = str(2:end-2);
str = fscanf(in,'figure_folder = %s',1);
figure_folder = str(2:end-2);
fclose(in);

%% Read config file

cfg = ReadConfig(Files);
folder_path   = [FE_folder cfg.output_folder];

%% Read data

filename_surf = getAllFiles(folder_path,'00_surface');
filename_moho = getAllFiles(folder_path,'01_surface');

data = load([folder_path 'physical_times.txt']);
t = data(:,2);

%% Hydrostatic equilibrium computation
% [fh,fval]=HydrostaticStateExact(r*1000,T,rho,0.1);

[fh,fval]=HydrostaticStateExact2l(...
    cfg.r_mean,...
    cfg.r_mean-cfg.depths_rho,...
    cfg.T,...
    cfg.rho(1),...
    cfg.rho(2),0.1, 0.1);

% outer shape
[a1,~,c1] = fr2abc(cfg.r_mean,fh(1),0);
% core
[a2,~,c2] = fr2abc(cfg.r_mean-cfg.depths_rho,fh(2),0);

fi = (-90:1:90);
lambda = (-180:1:180);
[fii,lambdai] = meshgrid(fi,lambda);

r1_ell = TriEllRadVec(fii/180*pi,lambdai/180*pi,a1,a1,c1,'rad');
r2_ell = TriEllRadVec(fii/180*pi,lambdai/180*pi,a2,a2,c2,'rad');

lmcosi_hydrostatic1 = xyz2plm(r1_ell',6);
lmcosi_hydrostatic2 = xyz2plm(r2_ell',6);

C20_1 = lmcosi_hydrostatic1(4,3);
C40_1 = lmcosi_hydrostatic1(11,3);
C60_1 = lmcosi_hydrostatic1(22,3);

C20_2 = lmcosi_hydrostatic2(4,3);
C40_2 = lmcosi_hydrostatic2(11,3);
C60_2 = lmcosi_hydrostatic2(22,3);


%% Figure for the movie

% draw and record all other frames

progressbar(0);

for i=1:numel(filename_surf)
    
    [~,name,~] = fileparts(filename_surf{i});
    
    if isempty(str2num(name(5:7)))
        ind(i) = str2num(name(5:6));
    else
        ind(i) = str2num(name(5:7));
    end
    
    lmcosi_limb = quad2plm(filename_surf{i},L);
    lmcosi_moho = quad2plm(filename_moho{i},L);
    
    lmcosi_limb(4,3)  =  lmcosi_limb(4,3) - C20_1;
    lmcosi_limb(11,3) = lmcosi_limb(11,3) - C40_1;
    lmcosi_limb(22,3) = lmcosi_limb(22,3) - C60_1;
    
    lmcosi_moho(4,3)  =  lmcosi_moho(4,3) - C20_1;
    lmcosi_moho(11,3) = lmcosi_moho(11,3) - C40_1;
    lmcosi_moho(22,3) = lmcosi_moho(22,3) - C60_1;
    
    C_only_topo(i) = lmcosi_limb((n_only+1)*n_only/2+1,3);
    C_only_moho(i) = lmcosi_moho((n_only+1)*n_only/2+1,3);
    
    tm_ratio(i) = C_only_topo/C_only_moho;

    progressbar(i/(numel(filename_surf)-1));
    
end

progressbar(1);

figure; hold on;
plot(t,abs(tm_ratio),'ok');
ylim([0 5]);

figure; hold on;
plot(t,-C_only_topo,'or');
plot(t,-C_only_moho,'ob');





