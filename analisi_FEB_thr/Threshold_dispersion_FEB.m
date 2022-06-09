% Threshold dispersion analysis for each channel and ASIC (on FEB)

%% Global optimisation
clear all; clc;

% global vars
main_folder = '/home/lucaghislotti/Downloads/analisi_FEB_thr/Analisi aggregate_FEB/';
global_folder = 'global_optimisation/';
asic_folder = 'asic_optimisation/';
output_data_folder = 'output_data/';
plots_folder = 'plots/';
boxplots_folder = 'boxplots/';
corrections_folder = 'corrections/';
globalthrdisp_folder = 'global thr dispersion/';

N = 500;
nonfun = 0;

%classA or full
opt_method = 'full';

if (strcmp(opt_method,'classA'))
    %pick A-class
    matA = readmatrix([main_folder, 'ASIC_Classification_A.dat']);
    matA = matA(4,[1:10]);
    range = matA;
    N = length(range);
    final_folder = 'classA opt/';
else
    range = 1:N-1;
    final_folder = 'global opt/';
end

% check for folder to exist
if (~exist(main_folder,'dir'))
    mkdir(main_folder);
end
% for each asic, each channel and each pt extract a,b (erf coeffs) from fitting
for iasic=range
    number = sprintf('%03i',iasic);
    if (exist(['FEB_Test/MODULE_',number,'/1/analysis_matlab/ThresholdScan/fitParameters.dat']) && exist(['FEB_Test/MODULE_',number,'/1/analysis_matlab/TransferFunction/Low_Energy_Gain.dat']))
        % threshold scan fit data
        data = importdata(['FEB_Test/MODULE_',number,'/1/analysis_matlab/ThresholdScan/fitParameters.dat']);
        temp = data.data(:,4:5);
        % tf low gain (high charge section) data
        data2 = importdata(['FEB_Test/MODULE_',number,'/1/analysis_matlab/TransferFunction/Low_Energy_Gain.dat']);
        for i=0:31
            % cols peaking time, rows fthr (get a,b)
            a_fit(iasic+1).ch(i+1).mat = reshape(temp(64*i+1:64*(i+1),1),8,8);
            b_fit(iasic+1).ch(i+1).mat = reshape(temp(64*i+1:64*(i+1),2),8,8);
            % peaking time (get a parameter from tf that signals whether
            % the channel works, i.e. gain)
            cubic_gain(iasic+1).mat(i+1, :) = data2.data(i*8+1:(i+1)*8, 15);
            % remove channels that dont work
            cubic_gain(iasic+1).mat(isnan(cubic_gain(iasic+1).mat) | cubic_gain(iasic+1).mat <= 0) = nan;
        end
    else
        % if some info do not exist, skip that asic (setting to nan doesnt
        % influence further analyses
        a_fit(iasic+1).ch(1).mat = nan;
        b_fit(iasic+1).ch(1).mat = nan;
		nonfun = nonfun + 1;
    end
end

disp('structures created!')


%% Mean over every ASIC, ch, tau, fthr
% put all data in a big matrix to compute mean without outliers and non
% functioning channels
global_mean = 0;
temp = zeros(N * 32 * 8 * 8, 1);
temp_cub = zeros(N * 32 * 8 * 8, 1);
for iasic=range
    if (~isnan(a_fit(iasic+1).ch(1).mat))
        for i = 0:31
            for j = 0:7
                for k = 0:7
                    temp(32*8*8*iasic + 8*8*i + 8*j + k + 1) = a_fit(iasic+1).ch(i+1).mat(k+1, j+1);
                    temp_cub(32*8*8*iasic + 8*8*i + 8*j + k + 1) = cubic_gain(iasic+1).mat(i+1, j+1);
                end
            end
        end
    end
end

global_mean = nanmean(rmoutliers(temp(~isnan(temp_cub))));
disp('global mean calculated!')


%% FTHR which gives least dispersion globally
for iasic=range
    best_fthr(iasic+1).mat = nan(32,8);
    if (~isnan(a_fit(iasic+1).ch(1).mat))
        for i=0:31
            % getting the idx (fthr+1)
            [~, idx] = min(abs(a_fit(iasic+1).ch(i+1).mat - global_mean),[],1);
            % rows channels, cols peaking time
            best_fthr(iasic+1).mat(i+1,:) = idx-1;
        end
    end
end
disp('global optimisation done!')


%% Standard deviation over selected FTHR, ASIC, ch, tau
global_std = 0;
for iasic=range
    temp = nan(32,8);
    if (~isnan(a_fit(iasic+1).ch(1).mat))
        for i=0:31
            current_fthr = best_fthr(iasic+1).mat(i+1,:)+1;
            for j=0:7
                temp(i+1,j+1) = a_fit(iasic+1).ch(i+1).mat(current_fthr(j+1),j+1);
            end
        end
        global_std = global_std + mean(std(temp));
    end
end
global_std = global_std/(N-nonfun);
disp('statistics for global optimisation obtained!')


%% Optimisation for each ASIC
% Mean over every ch, tau, fthr for each ASIC
asic_mean_v = zeros(N,1);
temp = zeros(32*8*8, 1);
temp_cub = zeros(32*8*8, 1);
for iasic=range
    if (~isnan(a_fit(iasic+1).ch(1).mat))
        for i=0:31
            for j=0:7
                for k=0:7
                    temp(8*8*i+8*j+k+1) = a_fit(iasic+1).ch(i+1).mat(k+1, j+1);
                    temp_cub(8*8*i + 8*j + k + 1) = cubic_gain(iasic+1).mat(i+1, j+1);
                end
            end
        end
        asic_mean_v(iasic+1) = nanmean(rmoutliers(temp(~isnan(temp_cub))));
    end
end
disp('per-asic mean calculated!')


%% FTHR which gives least dispersion for each ASIC
for iasic=range
    best_fthr2(iasic+1).mat = nan(32,8);
    if (~isnan(a_fit(iasic+1).ch(1).mat))
        for i=0:31
            % getting the idx (fthr+1)
            [~, idx] = min(abs(a_fit(iasic+1).ch(i+1).mat - asic_mean_v(iasic+1)),[],1);
            % rows channels, cols peaking time
            best_fthr2(iasic+1).mat(i+1,:) = idx-1;
        end
    end
end
disp('per-asic optimisation done!')


%% Standard deviation over selected FTHR, ASIC, ch, tau
asic_std = zeros(N,1);
for iasic=range
    temp = nan(32,8);
    if (~isnan(a_fit(iasic+1).ch(1).mat))
        for i=0:31
            current_fthr = best_fthr2(iasic+1).mat(i+1,:)+1;
            for j=0:7
                temp(i+1,j+1) = a_fit(iasic+1).ch(i+1).mat(current_fthr(j+1),j+1);
            end
        end
        asic_std(iasic) = mean(std(temp));
    end
end
disp('statistics for per-asic optimisation obtained!')


%% Export
if (~exist(main_folder,'dir'))
    mkdir(main_folder);
end
if (~exist([main_folder, global_folder],'dir'))
    mkdir([main_folder, global_folder]);
end
if (~exist([main_folder, asic_folder],'dir'))
    mkdir([main_folder, asic_folder]);
end
if (~exist([main_folder, global_folder, output_data_folder, final_folder],'dir'))
    mkdir([main_folder, global_folder, output_data_folder, final_folder]);
end
if (~exist([main_folder, asic_folder, output_data_folder],'dir'))
    mkdir([main_folder, asic_folder, output_data_folder]);
end
if (~exist([main_folder, global_folder, plots_folder],'dir'))
    mkdir([main_folder, global_folder, plots_folder]);
end
if (~exist([main_folder, global_folder, plots_folder, boxplots_folder],'dir'))
    mkdir([main_folder, global_folder, plots_folder, boxplots_folder, final_folder]);
end
if (~exist([main_folder, global_folder, plots_folder, corrections_folder],'dir'))
    mkdir([main_folder, global_folder, plots_folder, corrections_folder, final_folder]);
end
if (~exist([main_folder, asic_folder, plots_folder],'dir'))
    mkdir([main_folder, asic_folder, plots_folder]);
end
if (~exist([main_folder, asic_folder, plots_folder, boxplots_folder],'dir'))
    mkdir([main_folder, asic_folder, plots_folder, boxplots_folder]);
end
if (~exist([main_folder, asic_folder, plots_folder, corrections_folder],'dir'))
    mkdir([main_folder, asic_folder, plots_folder, corrections_folder]);
end
disp('folders created!')


%% Threshold dispersion aggregate analysis
% Conversion type
type = 'SS'; % FF,TT,SS

channels = 0:1:31;

for iasic=range
    
    % Global best fthr
    best_mat = best_fthr(iasic+1).mat;
    
    file_out = sprintf([main_folder, global_folder, output_data_folder, final_folder, 'ASIC%03i_best_fthr.dat'],iasic);
    
    fileID = fopen(file_out,'w');
    fprintf(fileID,'#---------------------------------------------------------\n');
    fprintf(fileID,'# Best fthr for global dispersion optimisation for ASIC%03i.\n',iasic);
    fprintf(fileID,'# Rows: #Channels\n');
    fprintf(fileID,'# Columns: Peaking time\n');
    fprintf(fileID,'#---------------------------------------------------------\n');
    s = repmat('%i\t',1,8);
    s(end+1:end+2) = '\n';
    fprintf(fileID,s,best_mat');
    fclose(fileID);
    
    % Plot for global and asic
    if (~isnan(a_fit(iasic+1).ch(1).mat))
        for j=0:7
            temp_gain = nan(32,1);
            for i=0:31                
                G_ch = cubic_gain(iasic+1).mat(i+1, j+1); % same for each one
                Thr_DAC = a_fit(iasic+1).ch(i+1).mat(4,j+1);
                
                % fixed fthr (for comparison w/ best fthr)
                Thr_011(iasic+1).mat(i+1, j+1) = conversion_DACthru_keV(G_ch, Thr_DAC, type);
                a_011(iasic+1).mat(i+1, j+1) = a_fit(iasic+1).ch(i+1).mat(4,j+1);
                
                % best fthr
                current_fthr = best_fthr(iasic+1).mat(i+1,j+1)+1;
                
                Thr_DAC = a_fit(iasic+1).ch(i+1).mat(current_fthr,j+1);
                
                Thr_best(iasic+1).mat(i+1, j+1) = conversion_DACthru_keV(G_ch, Thr_DAC, type);
                a_best(iasic+1).mat(i+1, j+1) = a_fit(iasic+1).ch(i+1).mat(current_fthr,j+1);

                % Minimised (asic)
                % best fthr
                current_fthr = best_fthr2(iasic+1).mat(i+1,j+1)+1;
                
                Thr_DAC = a_fit(iasic+1).ch(i+1).mat(current_fthr,j+1);
                
                Thr_best_single(iasic+1).mat(i+1, j+1) = conversion_DACthru_keV(G_ch, Thr_DAC, type);
                a_best_single(iasic+1).mat(i+1, j+1) = a_fit(iasic+1).ch(i+1).mat(current_fthr,j+1);
            end
        end
    else
        Thr_011(iasic+1).mat = nan(32, 8);
        Thr_best(iasic+1).mat = nan(32, 8);
        Thr_best_single(iasic+1).mat = nan(32, 8);
        a_011(iasic+1).mat = nan(32, 8);
        a_best(iasic+1).mat = nan(32, 8);
        a_best_single(iasic+1).mat = nan(32, 8);
    end
end
disp('fine threshold structures and file report created!')


%% Threshold dispersion globally before and after trimming per pt
clear acc_011 acc_opt

cmap_trans = [[51 66 255 256/4]/256;[255 51 252 256/4]/256];
cmap = [[51 66 255 256]/256;[255 51 252 256]/256];

acc_011 = 0;
acc_opt = 0;
for j = 0:7
    file_out = [main_folder, global_folder, plots_folder, globalthrdisp_folder, final_folder,  sprintf(['global_threshold_dispersion_tau%i_' opt_method '.pdf'],j)];
    
    % Global threshold dispersion
    idx_011 = 1;
    idx_opt = 1;

    if(strcmp(opt_method, 'classA') == 1)
        gbl_thr_range = matA;
    else
        gbl_thr_range = 1:N-1;
    end

    for iasic = gbl_thr_range
        if (iasic == 132)
            continue
        end
        
        Thr_011_no_out = Thr_011(iasic+1).mat(~isoutlier(Thr_011(iasic+1).mat(:, j+1)), j+1);
        Thr_no_out = Thr_best(iasic+1).mat(~isoutlier(Thr_best(iasic+1).mat(:, j+1)),j+1);
        
        %add sample to global arrays
        count_011 = length(Thr_011_no_out);
        count_opt = length(Thr_no_out);
        acc_011(idx_011:idx_011+count_011-1,1) = Thr_011_no_out';
        acc_opt(idx_opt:idx_opt+count_opt-1,1) = Thr_no_out';
        idx_011 = idx_011 + count_011;  
        idx_opt = idx_opt + count_opt;
    end
    
    % plot trivial (011) vs trimmed
    f = figure('visible','off');
    set(f, 'Position', [10 10 500 400])
    pd_011 = fitdist(acc_011,'Normal');
    pd_opt = fitdist(acc_opt,'Normal');
    hold on;

    %hist
    yyaxis left
    h_011_hist = histogram(acc_011,'BinWidth',1);
    h_011_hist.FaceColor = cmap_trans(1,1:3);
    h_011_hist.FaceAlpha = cmap_trans(1,4);
    h_opt_hist = histogram(acc_opt,'BinWidth',1);
    h_opt_hist.FaceColor = cmap_trans(2,1:3);
    h_opt_hist.FaceAlpha = cmap_trans(2,4);

    %fit
    yyaxis right
    step = 0.01;
    limits = get(gca, 'XLim');
    pd = fitdist(acc_011, 'Normal');    % fit the data with normal distribution
    pd2 = makedist('Normal', pd.mu, pd.sigma);    % make the pdf
    y = pdf(pd2, limits(1):step:limits(2));
    y = y/sum(y);
    plot(limits(1):step:limits(2),y,'LineWidth',0.8,'Color',cmap(1,:),'LineStyle','-');

    pd = fitdist(acc_opt, 'Kernel');    % fit the data with kernel distribution
    % no need to make an object
    y = pdf(pd, limits(1):step:limits(2));
    y = y/sum(y);
    plot(limits(1):step:limits(2),y,'LineWidth',0.8,'Color',cmap(2,:),'LineStyle','-');

    ax = gca;
    ax.YAxis(2).Visible = 'off';
    box on;
    axis tight;
    
    xlabel('Threshold [keV]');
    ylabel('Occurrencies');
    title(sprintf('Global threshold dispersion for peaking time #%i', j));
    notes = {['Occurrencies_{011}: ',sprintf('%i',length(acc_011(~isnan(acc_011))))], ...
             ['\mu_{011}: ',sprintf('%1.3f',pd_011.mu),' [keV]'], ...
             ['\sigma_{011}: ',sprintf('%1.3f',pd_011.sigma),' [keV]'], ...
             ['Occurrencies_{opt}: ',sprintf('%i',length(acc_opt(~isnan(acc_opt))))], ...
             ['\mu_{opt}: ',sprintf('%1.3f',pd_opt.mu),' [keV]'], ...
             ['\sigma_{opt}: ',sprintf('%1.3f',pd_opt.sigma),' [keV]']
             };
    annotation('textbox','Position',[0.6 0.7 0.2 0.2],'String',notes,'FitBoxToText','on','BackgroundColor','white','FontSize',8);
    exportgraphics(f,file_out,'ContentType','vector');
    file_fig_out = [main_folder, global_folder, plots_folder, globalthrdisp_folder, final_folder,  sprintf(['global_threshold_dispersion_tau%i_' opt_method '.fig'],j)];
    savefig(file_fig_out);
    close(f);
end
disp('global threshold distribution plot created!')


%% Mean over every ASIC, ch, tau, fthr (Thr_011 [keV])
% put all data in a big matrix to compute mean without outliers and non
% functioning channels
global_mean_011_kev = 0;
temp = NaN(length(range)*32*8, 1);

count = 1;
for iasic=range
    if (~isnan(Thr_011(iasic).mat))
        for i = 1:32
            for j = 1:8
                if(Thr_011(iasic).mat(i, j)~=0 && ~isnan(Thr_011(iasic).mat(i, j)))
                    temp(count) = Thr_011(iasic).mat(i, j);
                    count = count + 1;
                end
            end
        end
    end
end

temp = rmoutliers(temp);
global_mean_011_kev = nanmean(temp);
disp('global mean calculated! [thr_kev]')