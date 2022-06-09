%% Threshold dispersion globally (not optimized)
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

    if(strcmp(opt_method, 'classA')==1)
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

    %fit
    yyaxis right
    step = 0.01;
    limits = get(gca, 'XLim');
    pd = fitdist(acc_011, 'Normal');    % fit the data with normal distribution
    pd2 = makedist('Normal', pd.mu, pd.sigma);    % make the pdf
    y = pdf(pd2, limits(1):step:limits(2));
    y = y/sum(y);
    plot(limits(1):step:limits(2),y,'LineWidth',0.8,'Color',cmap(1,:),'LineStyle','-');

    ax = gca;
    ax.YAxis(2).Visible = 'off';
    box on;
    axis tight;
    
    xlabel('Threshold [keV]');
    ylabel('Occurrencies');
    title(sprintf('Global threshold dispersion for \\tau_{%i}', j));
    notes = {['Occurrencies_{011}: ',sprintf('%i',length(acc_011(~isnan(acc_011))))], ...
             ['\mu_{011}: ',sprintf('%1.3f',pd_011.mu),' [keV]'], ...
             ['\sigma_{011}: ',sprintf('%1.3f',pd_011.sigma),' [keV]']
             };
    annotation('textbox','Position',[0.6 0.7 0.2 0.2],'String',notes,'FitBoxToText','on','BackgroundColor','white','FontSize',8);
    exportgraphics(f,file_out,'ContentType','vector');
    close(f);
end
disp('global threshold distribution plot created!')