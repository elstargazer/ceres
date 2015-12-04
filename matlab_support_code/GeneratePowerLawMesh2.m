function GeneratePowerLawMesh2(Files,cfg, Nrand)


matlab_config_filename = Files.matlab_config_filename;
config_template_filename = Files.config_template_filename;
config_list_filename = Files.config_list_filename;

%% Input paramters

r_mean    = cfg.r_mean;
beta      = cfg.beta;
intercept = cfg.intercept;
layer_mat = cfg.mat_id;

L = 50;

nsq = 24;
nl  = [30 15];

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

for i=1:Nrand
    
    lmcosi_shape = PowerLawSH(r_mean,beta,intercept,L); 
    
    % make lmcosi_cmb = [0 0 r_mean-cfg.depths_rho 0] to make a spherical
    % core
    lmcosi_cmb = PowerLawSH(r_mean-cfg.depths_rho,beta-1,intercept,L);
  
    % make shape always oblate
    lmcosi_shape(4,3) = -abs(lmcosi_shape(4,3));

    % make core always oblate
    lmcosi_cmb(4,3) = -abs(lmcosi_cmb(4,3));
    
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


