ccc
fclose all;

GM = 62.68e9;
G  = 6.67e-11;
Rref = 470000;
M = GM/G;

runname = 'benchmark_g30'; %mod this line to reflect the set of runs being computed
runs_root_folder = '/Users/rogerfu/Dropbox/ceres_dropbox/ceres_public/';
%runlist_filename file must be in the runs_root_folder
runlist_filename = [runs_root_folder 'RunList_benchmark_g30.txt']; %mod this line for the runlist .txt file
in_runlist = fopen(runlist_filename,'r');
output_general_folder = [runs_root_folder 'surface_outputs/'];

L = 100;
 
config_filename = fgetl(in_runlist);

all_spectra_folder = [output_general_folder runname '/spectra_' runname ];
mkdir(all_spectra_folder);

output_number = 1;

% figure; hold on;

while (config_filename ~= -1)
    
    % read config file
    Files.config_template_filename = [runs_root_folder config_filename];
    cfg = ReadConfig(Files);
    
    %     filename_surf = getAllFiles(folder_path,'_surface');
    %     filename_mesh = getAllFiles(folder_path,'_mesh');
    %
    %     data = load([folder_path '/physical_times.txt']);
    %     t = data(:,2);
    
    %% Hydrostatic equilibrium computation
    
    M1 = 4/3*pi*(cfg.r_mean)^3*cfg.rho(1);
    M2 = 4/3*pi*(cfg.r_mean - cfg.depths_rho)^3*(cfg.rho(2)-cfg.rho(1));
    
    [fh,fval]=HydrostaticStateExact2l(...
        cfg.r_mean,...
        cfg.r_mean-cfg.depths_rho,...
        cfg.T,...
        cfg.rho(1),...
        cfg.rho(2),0.1, 0.1);
    
    [a1,c1]=f2axes(cfg.r_mean,fh(1));
    
    fi = (-90:1:90);
    lambda = (-180:1:180);
    [fii,lambdai] = meshgrid(fi,lambda);
    r_ell = TriEllRadVec(fii/180*pi,lambdai/180*pi,a1,a1,c1,'rad');
    
    lmcosi_hydrostatic1 = xyz2plm(r_ell',6);
    
    C20_1 = lmcosi_hydrostatic1(4,3);
    C40_1 = lmcosi_hydrostatic1(11,3);
    C60_1 = lmcosi_hydrostatic1(22,3);
    
    
    %% read time data
    
    phys_times_data = load([output_general_folder runname '/output_' num2string(output_number,'%iu') '/physical_times.txt']);
    t = phys_times_data(:,2);
    
    
    %% compute J2 and fp1
%     ellipse_fits_file = [output_general_folder runname '/output_' num2string(output_number,'%iu') '/ellipse_fits.txt'];
%     
%     in_ell = fopen(ellipse_fits_file,'r');
%     ell_data = textscan(in_ell,'%d a_%d = %d  c_%d = %d\n');
%     fclose(in_ell);
% 
%     t_ind = ell_data{1}+1;
%     a_ind = ell_data{2}+1;
%     
%     clear ell1 ell2
% 
%     ell1(:,1) = double(ell_data{3}(a_ind == 1));
%     ell1(:,2) = double(ell_data{5}(a_ind == 1));
%     
%     ell2(:,1) = double(ell_data{3}(a_ind == 2));
%     ell2(:,2) = double(ell_data{5}(a_ind == 2));
%     
%     fp1 = (ell1(:,1)-ell1(:,2))./ell1(:,1);
%     
%     J2=-((ell1(:,2).*ell1(:,2)-ell1(:,1).*ell1(:,1)).*M1+...
%         ((ell2(:,2).*ell2(:,2)-ell2(:,1).*ell2(:,1)).*M2))./(5*Rref*Rref.*M)/sqrt(5);
%     
%     C1 = 0.4*ell1(:,1).^2.*M1;
%     C2 = 0.4*ell2(:,1).^2.*M2;
%     C = C1 + C2;
%     
%     NMOI = C./(M.*Rref^2);
%     
%     omega_init = 2*pi/(cfg.T*3600);
%     
%     L_init = C(1)*(omega_init);
%     
%     omega = omega_init * C(1)./C;
%     Ti = 2*pi./omega/3600;
  
%    plot J2 and fp1
%     plot(t,fp1(1:end-1)./fp1(1),'r'); 
%     plot(t,J2(1:end-1)./J2(1),'b');
%     drawnow

    %% Compute and write spectral power density
    
    spectra_folder = [all_spectra_folder '/output_' num2str(output_number)];
    mkdir(spectra_folder);
    disp(['Writing to folder ' spectra_folder]);
    
    filename_surf = getAllFiles([output_general_folder runname '/output_' num2string(output_number,'%iu') '/'],'00_surface');
    
%     figure; hold on;
%     set(gca,'XScale','log');
%     set(gca,'YScale','log');
%     box on;
    
     
        for i=1:numel(filename_surf)
            
             try
            
            [path,name,ext] = fileparts(filename_surf{i});
            
            output_spectrum_filename = [output_general_folder runname '/output_' num2string(output_number,'%iu') '/' strrep(name,...
                '00_surface', 'spectrum.txt')];
            
            % get outer boundary spherical harmonic expansion
            lmcosi_limb = quad2plm(filename_surf{i},L);
            
            % get inner boundary spherical harmonic expansion
            %lmcosi_limb_inner = quad2plm(filename_surf_inner{i},L);
            
%             w = [M1/M M2/M];
    
%             lmcosi_grav_fe = WeightSumExpansion(w,{lmcosi_limb,lmcosi_limb_inner});
            
            % subtract hydrostatic signal from topography
            lmcosi_limb(4,3) = lmcosi_limb(4,3) - C20_1;
            lmcosi_limb(11,3) = lmcosi_limb(11,3) - C40_1;
            lmcosi_limb(22,3) = lmcosi_limb(22,3) - C60_1;
            
            % subtract hydrostatic signal from gravity
%             lmcosi_grav_fe(4,3) = lmcosi_grav_fe(4,3) - C20_1;
%             lmcosi_grav_fe(11,3) = lmcosi_grav_fe(11,3) - C40_1;
%             lmcosi_grav_fe(22,3) = lmcosi_grav_fe(22,3) - C60_1;
            
            if(i==1)
                initial_lmcosi = lmcosi_limb;
            end
            
            [sdl_limb,l_limb] = plm2spec(lmcosi_limb);
%             [sdl_grav,l_grav] = plm2spec(lmcosi_grav);
            
            % record spectra in output default folder
            in_spec = fopen(output_spectrum_filename, 'w');
            
            fprintf(in_spec, '%2.8E\n', t(i));
            fprintf(in_spec, '%4i %2.8E\n', [l_limb sdl_limb/1e6]');
            
            fclose(in_spec);
            
             %% record topography spectra in special folder
            output_spectrum_filename_topo = [spectra_folder '/' strrep(name,...
                'surface', 'spectrum_topo.txt')];
            
            in_spec = fopen(output_spectrum_filename_topo, 'w');
            fprintf(in_spec, '%2.8E\n', t(i));
            fprintf(in_spec, '%4i %2.8E\n', [l_limb sdl_limb/1e6]');
            fclose(in_spec);
            
            %% record gravity spectra in special folder          
%             output_spectrum_filename_grav = [spectra_folder '/' strrep(name,...
%                 'surface', 'spectrum_grav.txt')];
%             
%             in_spec = fopen(output_spectrum_filename_grav, 'w');
%             fprintf(in_spec, '%2.8E\n', t(i));
%             fprintf(in_spec, '%4i %2.8E\n', [l_grav sdl_grav/1e6]');
%             fclose(in_spec);

            catch
                disp('error');
            end

         end
    
    config_filename = fgetl(in_runlist);
    output_number = output_number + 1;
    
end

fclose(in_runlist);
% close all





