function GeneratePowerLawMesh2(Files,cfg, Nrand)


matlab_config_filename = Files.matlab_config_filename;
config_template_filename = Files.config_template_filename;
config_list_filename = Files.config_list_filename;

%% Input paramters

r_mean    = cfg.r_mean;
beta      = cfg.beta;
intercept = cfg.intercept;
layer_mat = cfg.mat_id;

L = 80;

nsq = 9*4;
nl  = [7*4 3*4];

%% plume

% plume_size = [100000];
% plume_r    = [0];
% plume_lat    = [20];
% plume_rho    = [500];
% plume_cell_mat = [2];
% [xp,~,zp] = sph2cart(0,plume_lat/180*pi,plume_r);

%% Read configuration file

in = fopen(matlab_config_filename);

str = fscanf(in,'spherical_mesh_filename = %s\n',1);
init_mesh_filename=str(2:end-2);
str = fscanf(in,'shape_folder = %s\n',1);
shape_folder = str(2:end-2);
str = fscanf(in,'shape_filename = %s\n',1);
shape_filename = str(2:end-2);
str = fscanf(in,'figure_folder = %s\n',1);
figure_folder = str(2:end-2);
str = fscanf(in,'figure_spetrum_filename = %s',1);
figure_spetrum_filename = str(2:end-2);

fclose(in);

%% plotting settings
FigureSettings

%% Initial parameters

cell_type = 'quad';
[path,name,ext] = fileparts(init_mesh_filename);

%% Generate random power law spectrum

% layer_mat(1) -> core
% layer_mat(0) -> outer shell


%% Compute hydrostatic shape

[fh,fval]=HydrostaticStateExact2l(...
    cfg.r_mean,...
    cfg.r_mean-cfg.depths_rho,...
    cfg.T,...
    cfg.rho(1),...
    cfg.rho(2),0.1, 0.1);

% outer shape
[a,~,c] = fr2abc(cfg.r_mean,fh(1),0);

fi = (-90:1:90);
lambda = (-180:1:180);
[fii,lambdai] = meshgrid(fi,lambda);

r_ell = TriEllRadVec(fii/180*pi,lambdai/180*pi,a,a,c,'rad');

lmcosi_hydrostatic1 = xyz2plm(r_ell',6);

C20_1 = lmcosi_hydrostatic1(4,3);
C40_1 = lmcosi_hydrostatic1(11,3);
C60_1 = lmcosi_hydrostatic1(22,3);

% core

[a,~,c] = fr2abc(cfg.r_mean-cfg.depths_rho,fh(2),0);
r_ell = TriEllRadVec(fii/180*pi,lambdai/180*pi,a,a,c,'rad');

lmcosi_cmb = xyz2plm(r_ell',6);

for i=1:Nrand
    
    % non hydrostatic part
    lmcosi_shape = PowerLawSH(r_mean,beta,intercept,L);
    
    % add hydrostatic part
    lmcosi_shape(4,3) = lmcosi_shape(4,3) + C20_1;
    lmcosi_shape(11,3) = lmcosi_shape(11,3) + C40_1;
    lmcosi_shape(22,3) = lmcosi_shape(22,3) + C60_1;
        
    meshStruct_def_quad = GenerateQuadLayerMesh(...
        lmcosi_cmb,lmcosi_shape,layer_mat,nsq,nl);
    figure; hold on;
    plot(meshStruct_def_quad.V(:,1),meshStruct_def_quad.V(:,2),'.');
    

%     for k=1:numel(plume_cell_mat)
%         % put a plume
%         plum_distance = sqrt( (meshStruct_def_quad.V(:,1)-xp(k)).^2 + ...
%             (meshStruct_def_quad.V(:,2)-zp(k)).^2 );
%         ind_plume = find(plum_distance<plume_size(k));
%         
%         [~,~,ib] = intersect(ind_plume,meshStruct_def_quad.E);
%         
%         
%         plot(meshStruct_def_quad.V(meshStruct_def_quad.E(ib),1),...
%             meshStruct_def_quad.V(meshStruct_def_quad.E(ib),2),'.r','MarkerSize',14);
%         meshStruct_def_quad.cell_mat(ib) = plume_cell_mat(k);
%         
%     end
    
    deformed_mesh_quad_filename = [path '/' name '_def_quad_' num2str(i) ext];
    
    FillConfigTemplate(config_template_filename,deformed_mesh_quad_filename,num2str(i))
    Write_ucd(meshStruct_def_quad,deformed_mesh_quad_filename,cell_type)
    
end


