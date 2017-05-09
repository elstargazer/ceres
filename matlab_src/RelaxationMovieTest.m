ccc

%% Image settings
fntsize = 20;
fntsize_sm = 10;
im_size=[0 0 20 20];
fig_folder='~/Dawn/Figures/';

ccj = {[1 0 0 ],...
    [0 1 0],...
    [0 0 1],...
    [0.6 0.6 0.6],...
    [0.6 0.6 0.6],...
    [0.6 0.6 0.6]};

%% Input parameters

L = 100;
Rref = 470000;

%% Files
matlab_config_filename   = '~/Dawn/FE/config/ConfigurationMatlab.cfg';
config_template_filename = '~/Dawn/FE/config/May5_g10_softint/May5_g10_softint_1.cfg';
runname                  = 'May5_g10_softint';
filename_real_spc        = '~/Dawn/FE/spectrum/PSD_SPG_HAMO_NH.txt';
n_only = 4;

%% Read configuration

movie_filename =[runname '.avi'];

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

filename_mesh = getAllFiles(folder_path,'_mesh');
filename_surf = getAllFiles(folder_path,'00_surface');
filename_moho = getAllFiles(folder_path,'01_surface');
filename_pl   = getAllFiles(folder_path,'_failurelocations00');
filename_ps   = getAllFiles(folder_path,'_principalstresses00');
filename_viscbase  = getAllFiles(folder_path,'_baseviscosities');
filename_viscreg   = getAllFiles(folder_path,'_viscositiesreg00');
filename_stress    = getAllFiles(folder_path,'_stresstensor00');
filename_flow      = getAllFiles(folder_path,'_flow00');

data = load([folder_path 'physical_times.txt']);
t = data(:,2);

% compute J2
% ell1 = data(:,4:5);
% ell2 = data(:,6:7);
% 
% fp1 = (ell(:,1)-ell(:,2))./ell(:,1);
% 
% J2=-((ell1(:,2).*ell1(:,2)-ell1(:,1).*ell1(:,1)).*M1+...
%     ((ell2(:,2).*ell2(:,2)-ell2(:,1).*ell2(:,1)).*M2))./(5*Rref*Rref.*M)/sqrt(5);

%% Load real spectrum

sdl_real = load(filename_real_spc);
l_real = (0:numel(sdl_real)-1)';

lambda_real=2*pi./l_real;
lambda_linear=lambda_real*(cfg.r_mean/1000);
k_real=1./lambda_linear;

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

ang = linspace(0,pi/2,100);

xell = a1*cos(ang);
zell = c1*sin(ang);

xell_cmb = a2*cos(ang);
zell_cmb = c2*sin(ang);

%% Figure for the movie

relax_pl = figure('Position',[1 1 1.2*1200 1.2*500],'Color','w');
pl_shape = subplot(1,2,2);
hold on;
set(gca, 'FontSize',fntsize);

xlabel('x [km]','FontSize',fntsize,'interpreter','latex');
ylabel('z [km]','FontSize',fntsize,'interpreter','latex');
box on;

xlim([0 600]);
ylim([0 600]);

% plot spectrum
pl_spectrum = subplot(1,2,1);
hold on; box on; grid on;
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'FontSize',fntsize);
xlabel('Frequency [cycles/km]','FontSize',fntsize,'interpreter','latex');
ylabel('Topography non-hydrostatic PSD [$\textrm{km}^{2}$]','FontSize',fntsize,...
    'interpreter','latex');

ylim([1e-6 1e2]);
xlim([3e-4 1e-1]);

lmcosi_limb = quad2plm(filename_surf{1},L);
lmcosi_moho = quad2plm(filename_moho{1},L);


lmcosi_limb(4,3)  =  lmcosi_limb(4,3) - C20_1;
lmcosi_limb(11,3) = lmcosi_limb(11,3) - C40_1;
lmcosi_limb(22,3) = lmcosi_limb(22,3) - C60_1;

lmcosi_moho(4,3)  =  lmcosi_moho(4,3) - C20_2;
lmcosi_moho(11,3) = lmcosi_moho(11,3) - C40_2;
lmcosi_moho(22,3) = lmcosi_moho(22,3) - C60_2;

% sdl_real(3) = sdl_real(3) - C20_1^2/5;
% sdl_real(5) = sdl_real(5) - C40_1^2/9;
% sdl_real(7) = sdl_real(7) - C60_1^2/13;

tm_ratio(1) = lmcosi_limb((n_only+1)*n_only/2+1,3)/...
              lmcosi_moho((n_only+1)*n_only/2+1,3);

[sdl_limb,l_limb] = plm2spec(lmcosi_limb);
[sdl_moho,l_moho] = plm2spec(lmcosi_moho);

lambda=2*pi./l_limb;
lambda_linear=lambda*(cfg.r_mean/1000);
kf=1./lambda_linear;

% real spec
h_real_spec = plot(k_real(3:end),sdl_real(3:end)/1e6,...
    '-o','MarkerSize',2,'Color','b','LineWidth',3);

% PrintWhite(relax_pl,[fig_folder 'relax_1.jpg']);

x_for_fit = log10(kf(3:2:end));
y_for_fit = log10(sdl_limb(3:2:end)/1e6);

p_init = polyfit(x_for_fit,y_for_fit,1);

h_fit_spec = plot(10.^[5; x_for_fit; -5],10.^(polyval(p_init,[5; x_for_fit; -5])),...
    '-o','MarkerSize',2,'Color','r','LineWidth',3);

plot_limb = plot(kf(1:2:end),sdl_limb(1:2:end)/1e6,...
    '-o','MarkerSize',4,'Color','k','LineWidth',1);

plot_moho = plot(kf(1:2:end),sdl_moho(1:2:end)/1e6,...
    '--*','MarkerSize',4,'Color','k','LineWidth',1);

legend([h_real_spec h_fit_spec plot_limb plot_moho],{'Observed','Power law fit','FE topo','FE moho'},'FontSize',fntsize);

% plot shape
subplot(pl_shape)
% fig_spec=figure('Color','w','Position',[1 1 1000 1000]);
% set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on; box on; axis equal;

xlabel('r [km]','FontSize',fntsize,'interpreter','latex');
ylabel('z [km]','FontSize',fntsize,...
    'interpreter','latex');

set(gca, 'FontSize',fntsize);
set(gca,'XTick',0:100:600);
set(gca,'YTick',0:100:600);

title(['t = ' num2str(t(1),'%6.2e') ' [y]'],'FontSize',fntsize,'interpreter','latex');

meshStruct = Read_ucd(filename_mesh{1});
V = meshStruct.V/1000;
E = meshStruct.E;
cell_mat = meshStruct.cell_mat;

% plot orig surface
surf_data = load([folder_path 'time001_00_surface.txt']);
h_orig_surf = plot(surf_data(:,1)/1000,surf_data(:,2)/1000,'k-','LineWidth',5);

xlim([0 600]);
ylim([0 600]);

%% draw and record first frame

delete(h_orig_surf)

for j=1:size(E,1)
    
    color = ccj{cell_mat(j)+1};
    p(j) = patch([V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)], ...
        [V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)],color);
    l(j) = line([V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)], ...
        [V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)],'Color','k');
end

plot(xell/1000,zell/1000,'y--','LineWidth',2);
plot(xell_cmb/1000,zell_cmb/1000,'y--','LineWidth',2);

% base viscocities
bv = load(filename_viscbase{2});
x = bv(:,1);
z = bv(:,2);

% plot faillure location
fl = load(filename_pl{1});
try
    failure_plot = plot(fl(:,1)/1000,fl(:,2)/1000,'xy','MarkerSize',10);
end
% plot stresses
s = load(filename_ps{2});
% p_stresses_plot = scatter(x/1000,z/1000,10,s(:,1)./s(:,2),'filled');
% caxis([0.2 5]);

% record first frame
v = VideoWriter(movie_filename);
open(v);

frame = getframe(gcf);
writeVideo(v,frame);

% draw and record all other frames
for i=2:1:numel(filename_mesh)-1
    
    meshStruct = Read_ucd(filename_mesh{i});
    V = meshStruct.V/1000;
    E = meshStruct.E;
    
    subplot(pl_shape);
    title(['t = ' num2str(t(i),'%6.2e') ' [y]'],'FontSize',fntsize,'interpreter','latex');
    
    for j=1:size(E,1)
        set(l(j),'XData',[V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)]);
        set(l(j),'YData',[V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)]);
        set(p(j),'XData',[V(E(j,1),1) V(E(j,2),1) V(E(j,3),1) V(E(j,4),1) V(E(j,1),1)]);
        set(p(j),'YData',[V(E(j,1),2) V(E(j,2),2) V(E(j,3),2) V(E(j,4),2) V(E(j,1),2)]);
    end
    
    try
        fl = load(filename_pl{i});
        set(failure_plot,'XData',fl(:,1)/1000,'YData',fl(:,2)/1000);
    end
    
    lmcosi_limb = quad2plm(filename_surf{i},L);
    lmcosi_moho = quad2plm(filename_moho{i},L);

    lmcosi_limb(4,3)  =  lmcosi_limb(4,3) - C20_1;
    lmcosi_limb(11,3) = lmcosi_limb(11,3) - C40_1;
    lmcosi_limb(22,3) = lmcosi_limb(22,3) - C60_1;
    
    lmcosi_moho(4,3)  =  lmcosi_moho(4,3) - C20_1;
    lmcosi_moho(11,3) = lmcosi_moho(11,3) - C40_1;
    lmcosi_moho(22,3) = lmcosi_moho(22,3) - C60_1;
    
    [sdl_limb,l_limb] = plm2spec(lmcosi_limb);
    [sdl_moho,l_moho] = plm2spec(lmcosi_moho);
    
    tm_ratio(i) = lmcosi_moho((n_only+1)*n_only/2+1,3)./...
                  lmcosi_limb((n_only+1)*n_only/2+1,3);
                  

    set(plot_limb,'YData',sdl_limb(1:2:end)/1e6);
    set(plot_moho,'YData',sdl_moho(1:2:end)/1e6);
     
    frame = getframe(gcf);
    writeVideo(v,frame);
    i/numel(filename_surf)
    
end

close all
close(v)

%% Plot evolution of moho and topo

figure; hold on;
plot(t(2:end),1./abs(tm_ratio(2:end)),'-ok');
plot(t(2:end),15*exp(-t(2:end)/0.1e12)+0.8,'-r');

ylim([0 5]);

s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,max(cdate)],...
               'Startpoint',[1 1]);
           
f = fittype('c+a*ext(-b/x','problem','n','options',s);

[c2,gof2] = fit(cdate,pop,f,'problem',2);

%% Fit evolution of moho and topo

% y_signal = abs(1./tm_ratio(2:end));
% x_time = t(2:end);
% 
% [parmhat, parmci] = expfit(x_time,y_signal,true(size(x_time)));
% 
% f = @(b,x) b(1).*exp(b(2).*x) + b(3);
% nrmrsd = @(b) norm(y_signal - f(b,x_time))  ;                        % Residual Norm Cost Function
% B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
% [B,rnrm] = fminsearch(nrmrsd, B0);  
% 
% x_plot = linspace(min(x_time), max(x_time));
% figure(1)
% plot(x_time, y_signal, 'pg')
% hold on
% plot(x_plot, f(B,x_plot), '-r')
% hold off
% grid
% xlabel('Time')
% ylabel('Signal')
% legend('Data', 'Fit', 'Location','NE')
% 
