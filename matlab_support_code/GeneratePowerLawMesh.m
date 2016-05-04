function GeneratePowerLawMesh(Files,r_mean, beta, intercept, Nrand)


matlab_config_filename = Files.matlab_config_filename;
config_template_filename = Files.config_template_filename;
config_list_filename = Files.config_list_filename;

%% Input paramters

L = 100;

% core axes
a_core = 410000;
c_core = 410000;

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

%% Read init sphere mesh

meshStruct = Read_ucd(init_mesh_filename);

x = meshStruct.V(:,1);
z = meshStruct.V(:,2);

%% Generate random power law spectrum

lmcosi_shape = CreateEmptylmcosi(L);
lmcosi_shape(1,3) = r_mean;
lmcosi_shape(2,3:4) = 0;
lmcosi_shape(3,3:4) = 0;

progressbar(0);

for i=1:Nrand
    
    for n=2:2:lmcosi_shape(end,1)
        lmcosi_shape((n+1)*n/2+1,3) = ((rand<0.5)*2-1)*sqrt((2*n+1)*...
            10^polyval([beta intercept],log10(n)));
    end
    
    % make it always oblate
    lmcosi_shape(4,3) = -abs(lmcosi_shape(4,3));
    
    %% Compute limb
    
    [lon,lat,r] = cart2sph(x,0,z);
    r_shape = plm2xyz(lmcosi_shape,lat*180/pi,lon*180/pi);
    
    r_new = r.*r_shape;
    [x_new, ~, z_new] = sph2cart(lon,lat,r_new);
    
    meshStruct_def.E = meshStruct.E;
    meshStruct_def.V = [x_new z_new zeros(size(x_new))];
    
    xv1 = meshStruct_def.V(meshStruct_def.E(:,1),1);
    xv2 = meshStruct_def.V(meshStruct_def.E(:,2),1);
    xv3 = meshStruct_def.V(meshStruct_def.E(:,3),1);
    xv4 = meshStruct_def.V(meshStruct_def.E(:,4),1);
    
    zv1 = meshStruct_def.V(meshStruct_def.E(:,1),2);
    zv2 = meshStruct_def.V(meshStruct_def.E(:,2),2);
    zv3 = meshStruct_def.V(meshStruct_def.E(:,3),2);
    zv4 = meshStruct_def.V(meshStruct_def.E(:,4),2);
    
    eps = 1e-5;
    ind_quad = find((xv1 > -eps) & (xv2 > -eps) & (xv3 > -eps) & (xv4 > -eps) & ...
        (zv1 > -eps) & (zv2 > -eps) & (zv3 > -eps) & (zv4 > -eps));
    
    meshStruct_def_quad.E = meshStruct_def.E(ind_quad,:);
    meshStruct_def_quad.V = meshStruct_def.V;
    meshStruct_def_quad.cell_mat = zeros(size(meshStruct_def_quad.E,1),1);
    meshStruct_def.cell_mat = zeros(size(meshStruct_def.E,1),1);
    
    r_cell_center = zeros(1,size(meshStruct_def_quad.E,1));   
    x_cell_center = zeros(1,size(meshStruct_def_quad.E,1));   
    y_cell_center = zeros(1,size(meshStruct_def_quad.E,1)); 
    z_cell_center = zeros(1,size(meshStruct_def_quad.E,1));
    
    for l = 1:size(meshStruct_def_quad.E,1)
        x_cell_center(l) = mean(meshStruct_def_quad.V(meshStruct_def_quad.E(l,:),1));
        y_cell_center(l) = mean(meshStruct_def_quad.V(meshStruct_def_quad.E(l,:),2));
        z_cell_center(l) = mean(meshStruct_def_quad.V(meshStruct_def_quad.E(l,:),3));
        
        r_cell_center(l) = sqrt(x_cell_center(l).*x_cell_center(l) + ...
            y_cell_center(l).*y_cell_center(l) + ...
            z_cell_center(l).*z_cell_center(l));
        
    end
    
%     plot(x_cell_center,y_cell_center,'.');
%     plot(meshStruct_def_quad.V(:,1),meshStruct_def_quad.V(:,2),'.');    
    ell_r = (x_cell_center.^2)./(a_core.^2) + (y_cell_center.^2)./(c_core.^2);
    core_condition = ell_r < 1.0;
    meshStruct_def_quad.cell_mat(core_condition)=1;
    
    deformed_mesh_filename = [path '/' name '_def_' num2str(i) ext];
    deformed_mesh_quad_filename = [path '/' name '_def_quad_' num2str(i) ext];
    
    FillConfigTemplate(config_template_filename,deformed_mesh_quad_filename,num2str(i))
    
    Write_ucd(meshStruct_def_quad,deformed_mesh_quad_filename,cell_type)
    Write_ucd(meshStruct_def,deformed_mesh_filename,cell_type)
    
    progressbar(i/Nrand);
    
end

progressbar(1);
