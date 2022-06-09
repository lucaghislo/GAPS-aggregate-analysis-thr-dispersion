%% New Optimization approach 
% Change to 'classA' in Threshold_dispersions_ASIC.m

Vth_m = pd_opt.mu;

% global thr bin width
Vth_gb = 4;

% ASICs means per tau
% init
asics_tau_means_kev = NaN(length(matA), 1);
asics_all_channels = [];
corrections_ASIC = zeros(length(matA), 2);
count_opt_arr = zeros(length(matA), 1);
corrections_ASIC = zeros(length(matA), 2);
asic_ID = NaN(length(matA), 1);

count = 1;
for iasic = matA
    if (iasic == 132)
        continue
    end
   
    % Thr_best for fthr data
    Thr_no_out = Thr_best(iasic+1).mat(~isoutlier(Thr_best(iasic+1).mat(:, 7)),7); 
   
    % outliers included 
    % Thr_no_out = Thr_best(iasic+1).mat(:,7);
    
    % add sample to global arrays
    % Thr_no_out = rmoutliers(Thr_no_out);
    count_opt_arr(count) = length(Thr_no_out);
    asic_ID(count) = iasic;
    asics_all_channels = [asics_all_channels; Thr_no_out];
    asics_tau_means_kev(count) = nanmean(nonzeros(Thr_no_out));

    classA_ASICs(count).ID = iasic;
    classA_ASICs(count).n_ch = length(Thr_no_out);
    classA_ASICs(count).mean = nanmean(nonzeros(Thr_no_out));
    classA_ASICs(count).channels = Thr_no_out;
    classA_ASICs(count).delta = 0;
    classA_ASICs(count).extraThrCh = 0;

    count = count +1;
end

Vth_ma_kev = asics_tau_means_kev;

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
count = 1;

for iasic=1:length(matA)
    stop = start + count_opt_arr(iasic)-1;

    disp([num2str(iasic), ': ', num2str(start), ':', num2str(stop)])
    
    acc_opt_temp(start:stop) = acc_opt_temp(start:stop) - Vth_gb * corrections_ASIC(iasic, 1);
    acc_opt_temp(start:stop) = acc_opt_temp(start:stop) + Vth_gb * corrections_ASIC(iasic, 2);

    start = stop + 1;
end

disp("correction done!")


%% plot channels after correction
display_plot(acc_opt_temp, 'kernel')


%% remove ASICs (1st method)
% obtain the best 250 ASICs

% order by the delta between the ASIC mean and the main mean
for i=1:length(classA_ASICs)
    classA_ASICs(i).delta = abs(classA_ASICs(i).mean - Vth_m)
end

classA_fields = fieldnames(classA_ASICs);
classA_cells = struct2cell(classA_ASICs);
classA_size = size(classA_cells);
classA_cells = reshape(classA_cells, classA_size(1), []);
classA_cells = classA_cells';

classA_cells = sortrows(classA_cells, 5, 'descend');

% 311 classA ASICs
display_plot(cell2mat(classA_cells(:, 4)), 'Kernel');

% 250 classA ASICs
display_plot(cell2mat(classA_cells(62:end, 4)), 'Kernel');


%% remove ASICs (2nd method)
% obstain the best 250 ASICs

classA_fields = fieldnames(classA_ASICs);
classA_cells = struct2cell(classA_ASICs);
classA_size = size(classA_cells);
classA_cells = reshape(classA_cells, classA_size(1), []);
classA_cells = classA_cells';

% build struct to match the ASIC ID and its channels
counter = 1;
for i=1:length(cell2mat(classA_cells(:, 1)))
    for j=1:classA_ASICs(i).n_ch-1
        classA_channels(counter).ID = classA_ASICs(i).ID;
        classA_channels(counter).chVal = classA_ASICs(i).channels(j);
        counter = counter + 1;
    end
end

for i=1:length(classA_channels)
    classA_channels(i).delta = abs(classA_channels(i).chVal - Vth_m);
end

classA_ch_fields = fieldnames(classA_channels);
classA_ch_cells = struct2cell(classA_channels);
classA_ch_size = size(classA_ch_cells);
classA_ch_cells = reshape(classA_ch_cells, classA_ch_size(1), []);
classA_ch_cells = classA_ch_cells';

classA_ch_cells = sortrows(classA_ch_cells, 3, 'descend');

for i=1:length(classA_ch_cells)
    if cell2mat(classA_ch_cells(i, 3))>Vth_gb
        ASIC_i_ID = cell2mat(classA_ch_cells(i, 1));
        classA_ASICs(find([classA_ASICs.ID] == ASIC_i_ID)).extraThrCh = classA_ASICs(find([classA_ASICs.ID] == ASIC_i_ID)).extraThrCh + 1;
    end
end

classA_fields = fieldnames(classA_ASICs);
classA_cells = struct2cell(classA_ASICs);
classA_size = size(classA_cells);
classA_cells = reshape(classA_cells, classA_size(1), []);
classA_cells = classA_cells';

classA_cells = sortrows(classA_cells, 6, 'descend');

% barplot to show the number of bad channels per ASIC
% figure
% bar([classA_ASICs.ID], sort([classA_ASICs.extraThrCh]))

bad_asics = cell2mat(classA_cells(:, 1));
bad_asics = bad_asics(1:61);

% remove ASICs with worst channels
for i=1:length(classA_ASICs)-1
    for count = 1:length(bad_asics)-1
        if classA_ASICs(i).ID == bad_asics(count)
            classA_ASICs(i).channels = NaN;
        end
    end
end

classA_fields = fieldnames(classA_ASICs);
classA_cells = struct2cell(classA_ASICs);
classA_size = size(classA_cells);
classA_cells = reshape(classA_cells, classA_size(1), []);
classA_cells = classA_cells';

% 250 classA ASICs
display_plot(cell2mat(classA_cells(:, 4)), 'Normal');


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