%% Data load
clear; clc;

load("saved_data/FEB_aggregated_ENC_tau6.mat")
load("saved_data/FEB_aggregated_gain_tau6.mat")
load("saved_data/FEB_aggregated_pedestal_tau6.mat")
load("saved_data/FEB_aggregated_thr_disp_tau6.mat")
load("saved_data/FEB_IDs.mat")


%% ENC vs Cubic Gain
figure
set(gcf,'position',[0,0,1920,1080])
sgtitle('ENC vs Cubic Gain for all FEBs channels')
for i = 1:32
    channel_gain = FEB_aggregated_gain_tau6(i, :)';
    channel_ENC = FEB_aggregated_ENC_tau6(i, :)';
    C = clusterdata([channel_gain, channel_ENC], 'maxclust', 10);
    subplot(4, 8, i)
    scatter(channel_gain, channel_ENC, 109, C, 'filled', 'hexagram')
    title(sprintf('Channel #%i', i))
    xlabel('Cubic Gain [ADC_{u}/DAC_{inju}]');
    ylabel('ENC [keV]');
end
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/ENC_vs_gain.svg')
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/ENC_vs_gain.fig')


%% ENC vs Cubic Gain single plot
ENC = FEB_aggregated_ENC_tau6';
ENC = ENC(:)';
ENC = ENC';

gain = FEB_aggregated_gain_tau6';
gain = gain(:)';
gain = gain';

C = clusterdata([gain, ENC], 'maxclust', 1);

count = 1;
classification = NaN(length(FEB_IDs)*32, 1);
for i = 1:length(FEB_IDs)
    for j=1:32
        classification(count, 1) = j;%FEB_IDs(i);
        count = count + 1;
    end
end

figure
scatter(gain, ENC, 150, classification, '.');
title('ENC vs Cubic Gain for all FEBs channels')
ylabel('ENC [keV]');
xlabel('Cubic Gain [ADC_{u}/DAC_{inju}]');
cb = colorbar;
cb.Ticks = linspace(0, 32, 33);
cb.TickLabels = num2cell(0:32) ;
casix([1,32])

saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/ENC_vs_gain_single_plot.svg')
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/ENC_vs_gain_single_plot.fig')


%% Threshold dispersion vs Cubic Gain
figure
set(gcf,'position',[0,0,1920,1080])
sgtitle('Threshold Dispersion vs Cubic Gain for all FEBs channels')
for i = 1:32
    channel_gain = FEB_aggregated_gain_tau6(i, [1:34 36:end])';
    channel_thr_disp = FEB_aggregated_thr_disp_tau6(i, :)';
    C = clusterdata([channel_gain, channel_thr_disp], 'maxclust', 10);
    subplot(4, 8, i)
    scatter(channel_gain, channel_thr_disp, 109, C, 'filled', 'hexagram')
    title(sprintf('Channel #%i', i))
    xlabel('Cubic Gain [ADC_{u}/DAC_{inju}]');
    ylabel('Threshold Dispersion [keV]');
end
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/thr-disp_vs_gain.svg')
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/thr-disp_vs_gain.fig')


%% Pedestal vs Cubic Gain
figure
set(gcf,'position',[0,0,1920,1080])
sgtitle('Pedestal vs Cubic Gain for all FEBs channels')
for i = 1:32
    channel_gain = FEB_aggregated_gain_tau6(i, :)';
    channel_pedestal = FEB_aggregated_pedestal_tau6(i, :)';
    C = clusterdata([channel_gain, channel_pedestal], 'maxclust', 10);
    subplot(4, 8, i)
    scatter(channel_gain, channel_pedestal, 109, C, 'filled', 'hexagram')
    title(sprintf('Channel #%i', i))
    xlabel('Cubic Gain [ADC_{u}/DAC_{inju}]');
    ylabel('Pedestal [DAC_{thr} units]');
end
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/pedestal_vs_gain.svg')
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/pedestal_vs_gain.fig')


%% ENC vs Threshold Dispersion
figure
set(gcf,'position',[0,0,1920,1080])
sgtitle('ENC vs Threshold Dispersion for all FEBs channels')
for i = 1:32
    channel_ENC = FEB_aggregated_ENC_tau6(i, [1:34 36:end])';
    channel_thr_disp = FEB_aggregated_thr_disp_tau6(i, :)';
    C = clusterdata([channel_thr_disp, channel_ENC], 'maxclust', 10);
    subplot(4, 8, i)
    scatter(channel_thr_disp, channel_ENC, 109, C, 'filled', 'hexagram')
    title(sprintf('Channel #%i', i))
    xlabel('Threshold Dispersion [keV]');
    ylabel('ENC [keV]');
end
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/ENC_vs_thr_disp.svg')
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/ENC_vs_thr_disp.fig')


%% correlation coefficients between ENC and Cubic Gain
rho_array = NaN(32, 1);

for i=1:32
    channel_ENC = FEB_aggregated_ENC_tau6(i, :)';
    channel_gain = FEB_aggregated_gain_tau6(i, :)';
    rho = corr(channel_gain, channel_ENC);
    rho_array(i) = rho;
    disp(rho)
end

figure
plot(rho_array);
xlim([0 33])
title('Correlation between ENC and Cubic Gain')
xlabel('#Channel')
ylabel('Correlation Coefficient')
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/ENC_vs_gain_corr_plot.svg')
saveas(gcf, '/home/lucaghislotti/Downloads/analisi_FEB_thr/plots/FEB/correlation_study/ENC_vs_gain_corr_plot.fig')