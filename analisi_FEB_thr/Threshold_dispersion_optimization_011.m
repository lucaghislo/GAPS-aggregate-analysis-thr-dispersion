%% New Optimization approach 

Vth_m = pd_opt.mu;

% global thr bin width
Vth_gb = 4;

% ASICs means per tau
% init
asics_tau_means_kev = NaN(length(gbl_thr_range), 1);
asics_all_channels = [];
corrections_ASIC = zeros(length(gbl_thr_range), 2);
count_opt_arr = zeros(length(gbl_thr_range), 1);
corrections_ASIC = zeros(N, 2);

for iasic = gbl_thr_range
    if (iasic == 132)
        continue
    end
    
    % Thr_no_out = Thr_best(iasic+1).mat(~isoutlier(Thr_011(iasic+1).mat(:, 7)),7);
    Thr_no_out = Thr_best(iasic+1).mat(:,7);
    
    % remove outliers (channels)
    % Thr_no_out = rmoutliers(Thr_no_out);
    
    %add sample to global arrays
    count_opt_arr(iasic) = length(Thr_no_out);
    asics_all_channels = [asics_all_channels;Thr_no_out];
    asics_tau_means_kev(iasic) = nanmean(Thr_no_out);

end

Vth_ma_kev = asics_tau_means_kev;
%Vth_ma_kev = rmoutliers(Vth_ma_kev);

% global thr mean
Vth_m = nanmean(Vth_ma_kev)

display_plot(Vth_ma_kev, 'Normal')


%% differences
D_kev = Vth_ma_kev - Vth_m;
flag_counter = 0;
changed_1 = true;
changed_2 = true;

while (changed_1 || changed_2)
    disp('new iteration')

    changed_1 = false;
    changed_2 = false;

    for i = 1:length(D_kev)
        if D_kev(i) > Vth_gb
            disp(['ASIC ', num2str(i), ' corrected!']);
            disp(['before: ', num2str(Vth_ma_kev(i))])

            Vth_ma_kev(i) = Vth_ma_kev(i) - Vth_gb;
            corrections_ASIC(i, 1) = corrections_ASIC(i, 1) + 1;

            disp([' after: ', num2str(Vth_ma_kev(i))])
            disp(' ')

            changed_1 = true;
        elseif D_kev(i) < -Vth_gb
            disp(['ASIC ', num2str(i), ' corrected!']);
            disp(['before: ', num2str(Vth_ma_kev(i))])

            Vth_ma_kev(i) = Vth_ma_kev(i) + Vth_gb; 
            corrections_ASIC(i, 2) = corrections_ASIC(i, 2) + 1;

            changed_2 = true;

            disp([' after: ', num2str(Vth_ma_kev(i))])
            disp(' ')
        end
    end

    flag_counter = flag_counter + 1;
    D_kev = Vth_ma_kev - Vth_m;
end

disp(['iterazioni: ', num2str(flag_counter-1)])
display_plot(Vth_ma_kev, 'Normal')


%% plot channels before correction
display_plot(asics_all_channels, 'Kernel')


%% implement correction on channels

acc_opt_temp = asics_all_channels;
start = 1;

for iasic=gbl_thr_range
    stop = start + count_opt_arr(iasic)-1;

    disp([num2str(iasic), ': ', num2str(start), ':', num2str(stop)])
    
    acc_opt_temp(start:stop) = acc_opt_temp(start:stop) - Vth_gb * corrections_ASIC(iasic, 1);
    acc_opt_temp(start:stop) = acc_opt_temp(start:stop) + Vth_gb * corrections_ASIC(iasic, 2);

    start = stop + 1;
end

disp("correction done!")


%% plot channels after correction
display_plot(acc_opt_temp, 'kernel')


%% plot function
function display_plot(dataset, distr)
    % colors
    cmap_trans = [[51 66 255 256/4]/256;[255 51 252 256/4]/256];
    cmap = [[51 66 255 256]/256;[255 51 252 256]/256];

    % plot
    f = figure('visible','on');
    pd1 = fitdist(dataset, 'Normal');
    hold on;
    
    %hist
    yyaxis left
    h = histogram(dataset,'BinWidth',1);
    h.FaceColor = cmap_trans(1,1:3);
    h.FaceAlpha = cmap_trans(1,4);
    
    %fit
    yyaxis right
    step = 0.01;
    limits = get(gca, 'XLim');

    if(strcmp(distr, 'Normal') == 0)
        pd = fitdist(dataset, 'Kernel');    % fit the data with kernel distribution
        % no need to make an object
        y = pdf(pd, limits(1):step:limits(2));
        y = y/sum(y);
        plot(limits(1):step:limits(2),y,'LineWidth',1,'Color',cmap(2,:),'LineStyle','-');
    else
        pd = fitdist(dataset, 'Normal');  
        pd2 = makedist('Normal', pd.mu, pd.sigma);  
        y = pdf(pd2, limits(1):step:limits(2));
        y = y/sum(y);
        plot(limits(1):step:limits(2),y,'LineWidth',1,'Color',cmap(1,:),'LineStyle','-');
    end

    ax = gca;
    ax.YAxis(2).Visible = 'off';
    box off;
    axis tight;

    xlabel('Threshold [keV]');
    ylabel('Occurrencies');
    title(sprintf('Global threshold dispersion for \\tau_{6}', j));
    notes = {['Occurrencies: ',sprintf('%i',length(dataset(~isnan(dataset))))], ...
             ['\mu: ',sprintf('%1.3f',pd1.mu),' [keV]'], ...
             ['\sigma: ',sprintf('%1.3f',pd1.sigma),' [keV]'], ...
             };
    annotation('textbox','Position',[0.68 0.7 0.2 0.2],'String',notes,'FitBoxToText','on','BackgroundColor','white','FontSize',8);
end