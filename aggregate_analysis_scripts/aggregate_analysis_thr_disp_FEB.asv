%% Threshold dispersion for each channel per ASIC plotted as heatmap

% clear; clc;

N = 300;

main_folder = 'aggregated_thr_disp/';
if (~exist(main_folder,'dir'))
    mkdir(main_folder);
end

thr_map = nan(32,N);
class_A_idxs = nan(8,N);
class_B_idxs = nan(8,N);

for j=0:7
    for i=1:N-1
        data = Thr_011; % raw data
        single_FEB = data(i+1);
        thr_map(:,i+1) = single_FEB.mat(:, j+1);
    end

    thr_map(isnan(thr_map) | thr_map <= 0) = nan;

    % remove invalid FEBs
    plot_thr_map = abs(thr_map(:,~isnan(thr_map(1,:)))-global_mean_011_kev);
    asic_n = 0:1:N-1;
    plot_asic_n = asic_n(~isnan(thr_map(1,:)));

    figure('units','normalized','outerposition',[0 0 .9 .9], 'Visible','off');
    subaxis(1,1,1,'Padding', 0.04, 'Margin', 0.01);
    heatmap(plot_thr_map,'XDisplayLabels',plot_asic_n,'YDisplayLabels',0:1:31);
    colormap('jet');
    xlabel('#FEB');
    ylabel('#Channel');
    caxis([0,20])
    title(sprintf('Threshold dispersion for each channel in each FEB [keV]. Peaking time = tau\\_%i',j));
    saveas(gcf,sprintf([main_folder, 'thr_disp_tau%i.svg'],j));
end