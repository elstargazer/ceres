function GeneratePowerLawMesh2(Files,r_mean, beta, intercept, Nrand)


matlab_config_filename = Files.matlab_config_filename;
config_template_filename = Files.config_template_filename;
config_list_filename = Files.config_list_filename;

%% Input paramters

L = 50;

% core axes
h = 200000;

nsq = 15;
nl  = [10 5];

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

lmcosi_shape = CreateEmptylmcosi(L);
lmcosi_shape(1,3) = r_mean;
lmcosi_shape(2,3:4) = 0;
lmcosi_shape(3,3:4) = 0;

lmcosi_cmb = CreateEmptylmcosi(L);
lmcosi_cmb(1,3) = r_mean-h;
lmcosi_cmb(2,3:4) = 0;
lmcosi_cmb(3,3:4) = 0;

layer_mat = [0 1];

for i=1:Nrand
    
    for n=2:2:lmcosi_shape(end,1)
        lmcosi_shape((n+1)*n/2+1,3) = ((rand<0.5)*2-1)*sqrt((2*n+1)*...
            10^polyval([beta intercept],log10(n)));
    end
    
    for n=2:2:lmcosi_cmb(end,1)
        lmcosi_cmb((n+1)*n/2+1,3) = ((rand<0.5)*2-1)*sqrt((2*n+1)*...
            10^polyval([beta-1 intercept],log10(n)));
    end
    
    % make shape it always oblate
    lmcosi_shape(4,3) = -abs(lmcosi_shape(4,3));
    
    meshStruct_def_quad = GenerateQuadLayerMesh(...
        lmcosi_cmb,lmcosi_shape,layer_mat,nsq,nl);
    
    deformed_mesh_quad_filename = [path '/' name '_def_quad_' num2str(i) ext];
    
    FillConfigTemplate(config_template_filename,deformed_mesh_quad_filename,num2str(i))   
    Write_ucd(meshStruct_def_quad,deformed_mesh_quad_filename,cell_type)
    
end


