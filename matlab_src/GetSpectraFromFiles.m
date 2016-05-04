% ccc
SECSINYEAR = 365.2422*86400;

%% plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Figures/';

%%

spectra_folder = '~/Dawn/FE/output/run102_nosign/spectra/';
key_for_spectra = 'spectrum';

all_content = dir(spectra_folder);
dirFlags = [all_content.isdir];
subFolders = all_content(dirFlags);

% figure; hold on;
% set(gca,'XScale','log');
% set(gca,'YScale','log');

t_all = [];
sdl_all = [];

for i=3:numel(subFolders);
    
    pathname = [spectra_folder subFolders(i).name '/'];
    filename_quad = getAllFiles(pathname,key_for_spectra);
    
    for j=1:numel(filename_quad)
        
        in = fopen(filename_quad{j},'r');
        
        t = textscan(in,'%f\n',1);
        
        t_all = [t_all t{1}];
        
        data = textscan(in,'%d %f\n');
        l = data{1};
        sdl = data{2};
        sdl_all = [sdl_all sdl];
        
        %         plot(l(1:2:end),sdl(1:2:end),'-k');
        fclose(in);
    end
end

%% Plot coefficient for a specific degree


% l_to_plot = 4;
%
% ind = find(l==l_to_plot);
% figure; hold on;
% % set(gca,'XScale','log');
% set(gca,'YScale','log');
% plot(t_all,sdl_all(ind,:),'.k')


%% find relaxation times
figure; hold on;
set(gca,'YScale','log');

% degrees to fit
l_for_fit = 2:2:100;
for i = 1:numel(l_for_fit)
    
    ind = find(l==l_for_fit(i));

    %fit 1
    %     [p,S] = polyfit(t_all,log(sdl_all(ind,:)),1);
    %     tau(i) = -1/p(1);
     
    % fit 2
    %     fitOpt = fitoptions('Normalize','off');
    %     [fit_obj,gof,output] = fit(t_all',log(sdl_all(ind,:))','poly1',fitOpt);
    
    max_sdl = max((sdl_all(ind,:)));
    cond = sdl_all(ind,:)/max_sdl > 0.1;
    
    t_select = t_all(cond);
    sdl_select = sdl_all(ind,cond);
    
    [fit_obj,gof,output] = fit(t_select',log(sdl_select)','poly1');
    
    coeffs = coeffvalues(fit_obj);
    coeffs_confint = confint(fit_obj);
    a = coeffs(1);
    b = coeffs(2);
    
    tau(i) = -1/a;
    tau_down(i) = -1/coeffs_confint(1,1);
    tau_up(i)   = -1/coeffs_confint(2,1);
    
    plot(t_select,sdl_select,'.k');
    plot(t_select,exp(polyval([a b],t_select)),'.r');
    
    waitforbuttonpress();
    
end

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;
set(gca,'YScale','log');
set(gca,'XScale','log');

color = 'y';

plot(l_for_fit,tau/SECSINYEAR,'-o','MarkerSize',5,'LineWidth',3,'Color',color);
plot(l_for_fit,tau_down/SECSINYEAR,'-','MarkerSize',3,'LineWidth',1,'Color',color);
plot(l_for_fit,tau_up/SECSINYEAR,'-','MarkerSize',3,'LineWidth',1,'Color',color);

xlabel('SH degree','FontSize',fntsize);
ylabel('Relaxation time [y]','FontSize',fntsize);











