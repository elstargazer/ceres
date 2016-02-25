ccc
SECSINYEAR = 365.2422*86400;

%% plotting settings
fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Figures/';

%%

spectra_folder = '~/Dawn/FE/output/run102/spectra/';

all_content = dir(spectra_folder);
dirFlags = [all_content.isdir];
subFolders = all_content(dirFlags);

figure; hold on;
set(gca,'XScale','log');
set(gca,'YScale','log');

t_all = [];
sdl_all = [];

for i=3:numel(subFolders);
    
    pathname = [spectra_folder subFolders(i).name '/'];
    filename_quad = getAllFiles(pathname,'spectrum');
    
    for j=1:numel(filename_quad)
        
        in = fopen(filename_quad{j},'r');
        
        t = textscan(in,'%f\n',1);
        
        t_all = [t_all t{1}];
        
        data = textscan(in,'%d %f\n');
        l = data{1};
        sdl = data{2};
        
        sdl_all = [sdl_all sdl];
        
        plot(l(1:2:end),sdl(1:2:end),'-k');
        fclose(in);
        
    end
end

%%


l_to_plot = 20;

ind = find(l==l_to_plot);
figure; hold on;
% set(gca,'XScale','log');
set(gca,'YScale','log');
plot(t_all,sdl_all(ind,:),'.k')


%% find relaxation times


figure;

% degrees to fit
l_for_fit = 2:2:100;
for i = 1:numel(l_for_fit)
    
    %fit 1
    ind = find(l==l_for_fit(i));
    [p,S] = polyfit(t_all,log(sdl_all(ind,:)),1);
    tau(i) = -1/p(1);
    
    % fit 2
    %     fitOpt = fitoptions('Normalize','off');
    %     [fit_obj,gof,output] = fit(t_all',log(sdl_all(ind,:))','poly1',fitOpt);
    
    max_sdl = max((sdl_all(ind,:)));
    cond = sdl_all(ind,:)/max_sdl > 0.1;
    
    t_select = t_all(cond);
    sdl_select = sdl_all(ind,cond);
    
    [fit_obj,gof,output] = fit(t_select',sdl_select','exp1');
    
    coeffs = coeffvalues(fit_obj);
    coeffs_confint = confint(fit_obj);
    a = coeffs(1);
    b = coeffs(2);
    
    tau2(i) = -1/b;
    tau2_down(i) = -1/coeffs_confint(1,2);
    tau2_up(i)   = -1/coeffs_confint(2,2);
       
%     plot(fit_obj,t_select,sdl_select)   
%     waitforbuttonpress();
    
end

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

plot(l_for_fit,tau/SECSINYEAR,'-ok','MarkerSize',3,'LineWidth',3)
% plot(l_for_fit,tau2/SECSINYEAR,'-or','MarkerSize',3,'LineWidth',3)
% errorbar(l_for_fit,tau2/SECSINYEAR,...
%     tau2_down/SECSINYEAR,tau2_up/SECSINYEAR);

xlabel('Degree','FontSize',fntsize);
ylabel('Relaxation time [y]','FontSize',fntsize);











