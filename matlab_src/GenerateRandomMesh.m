ccc;

%% plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Figures/';

%% Initial parameters

r_sph = 470000;
cell_type = 'quad';

init_mesh_filename = '../FE/meshes/mesh_sph.inp';
[path,name,ext] = fileparts(init_mesh_filename);

shape_folder='/Users/antonermakov/Dawn/CeresShapeModel/SPC/CERES_SURVEY_150716_GRAVITY_SPC/';
shape_filename='SHAPE_SPC150716_256.bds';

filename_quad = '../CeresFE/FE/outer_vertices00.txt';

Nrand = 10;

%% Read init sphere mesh

meshStruct = Read_ucd(init_mesh_filename);

x = meshStruct.V(:,1);
z = meshStruct.V(:,2);

figure; hold on;
plot(x,z,'.');
axis equal;

real_shape_full_filename = [shape_folder shape_filename];

step = 1;
Npts = 100000;
[x_grid,y_grid,z_grid]=ReadSPC(real_shape_full_filename,step,'grid');

x_grid = x_grid*1000;
y_grid = y_grid*1000;
z_grid = z_grid*1000;

r_grid=sqrt(x_grid.^2+y_grid.^2+z_grid.^2);

figure; hold on;
surf(x_grid,y_grid,z_grid,r_grid); shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

lmcosi_real_shape = xyz2plm(flipud(r_grid'));
[r_sh, lon_sh, lat_sh] = plm2xyz(lmcosi_real_shape);
[lon_sh,lat_sh] = meshgrid(lon_sh,lat_sh);
[x_sh,y_sh,z_sh] = sph2cart(lon_sh/180*pi,lat_sh/180*pi,r_sh);

figure; hold on;
surf(x_sh,y_sh,z_sh,r_sh); shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

%% Ceres spectrum

[sdl,l] = plm2spec(lmcosi_real_shape);

fig_spec=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;
set(gca,'XScale','log');
set(gca,'YScale','log');

pl_real = plot(l,sdl,'-k');

xlabel('Degree','FontSize',fntsize);
ylabel('Power [km^{2}]','FontSize',fntsize);

gi = ginput(2);

minn=round(gi(1,1));
maxn=round(gi(2,1));

cond = ((l>minn) & (l<maxn));

l_c   = l(cond);
sdl_c = sdl(cond);

p = polyfit(log10(l_c),log10(sdl_c),1);
ipt1 = polyval(p,log10(1));

%% Generate random spectrum

L = 180;

lmcosi_shape = CreateEmptylmcosi(L);
lmcosi_shape(1,3) = lmcosi_real_shape(1,3);
lmcosi_shape(2,3:4) = 0;
lmcosi_shape(3,3:4) = 0;

for i=1:Nrand
    
    for n=2:2:lmcosi_shape(end,1)
        %     lmcosi_shape((n+1)*n/2+1,3) = ((rand<0.5)*2-1)*sqrt((2*n+1)*sdl_r(n+1));
        lmcosi_shape((n+1)*n/2+1,3) = ((rand<0.5)*2-1)*sqrt((2*n+1)*10^polyval(p,log10(n)));
    end
    
    % make it always oblate
    lmcosi_shape(4,3) = -abs(lmcosi_shape(4,3));
    
    %     [sdl_r,l_r] = plm2spec(lmcosi_shape);
    %
    %     figure(fig_spec);
    %     pl_fit = plot(l_r(1:2:end),sdl_r(1:2:end),'-og','MarkerSize',2);
    
    %% Compute limb
    
    [lon,lat,r] = cart2sph(x,0,z);
    r_shape = plm2xyz(lmcosi_shape,lat*180/pi,lon*180/pi);
    
    r_new = r.*r_shape;
    [x_new, y_new, z_new] = sph2cart(lon,lat,r_new);
    
    %     figure; hold on;
    %     plot(x_new,z_new,'.');
    
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
    
    deformed_mesh_filename = [path '/' name '_def_' num2str(i) ext];
    deformed_mesh_quad_filename = [path '/' name '_def_quad_ ' num2str(i) ext];
    
    Write_ucd(meshStruct_def_quad,deformed_mesh_quad_filename,cell_type)
    Write_ucd(meshStruct_def,deformed_mesh_filename,cell_type)
    
end

%% limb spectrum
L = 30;

lmcosi_limb = quad2plm(filename_quad,L);

[sdl_limb,l_limb] = plm2spec(lmcosi_limb);
figure(fig_spec);
pl_relax = plot(l_limb(1:2:end),sdl_limb(1:2:end)...
    ,'-ro','MarkerSize',2);
legend([pl_real pl_fit pl_relax],...
    {'Real Ceres','Fit','Relaxed'},'FontSize',fntsize_sm);

PrintWhite(fig_spec,[fig_folder 'Fig_RelaxedCeresSpectrum.jpg']);

[r_limb_sh, lon_limb_sh, lat_limb_sh] = ...
    plm2xyz(lmcosi_limb);

[lon_limb_sh,lat_limb_sh] = meshgrid(lon_limb_sh,lat_limb_sh);
[x_limb_sh,y_limb_sh,z_limb_sh] = ...
    sph2cart(lon_limb_sh/180*pi,lat_limb_sh/180*pi,r_limb_sh);

figure; hold on;
surf(x_limb_sh,y_limb_sh,z_limb_sh,r_limb_sh); shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

