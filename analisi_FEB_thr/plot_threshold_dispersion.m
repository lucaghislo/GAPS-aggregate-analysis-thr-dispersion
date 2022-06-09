index = 0;
counter = 0;

for num_folder = 1:500
    folder = strcat('/home/lucaghislotti/Downloads/analisi_FEB_thr/ASIC_Test/ASIC_',num2str(num_folder, '%03d'));
    if (exist(folder, 'dir'))
        counter = counter + 1;
    end
end

thr_disp = zeros(counter, 8);

for num_folder = 1:500
    
    folder = strcat('/home/lucaghislotti/Downloads/analisi_FEB_thr/ASIC_Test/ASIC_',num2str(num_folder, '%03d'));
    
    if (exist(folder, 'dir'))
        file_thr = strcat(folder, '/1/analysis_matlab/ThresholdScan/Threshold_dispersion.dat');
        data = importdata(file_thr);
        index = index + 1;

        for t = 1:8
            thr_disp(index, t) = data.data(t, 4);
        end

    end

end

 
for i = 1:8
    figure
    histfit(thr_disp(:, i), 10);
    norm = fitdist(thr_disp(:, i), 'Normal');
    ylim([0, 3])
    title({['Threshold dispersion for \tau = ', num2str(i-1)], ['\mu = ',  num2str(round(norm.mu, 3)),' keV, \sigma = ', num2str(round(norm.sigma, 3)), ' keV']})
    hold on
    ax = gca;
    exportgraphics(ax,['plots/thrdisp_',num2str(i-1),'.pdf'],'Resolution',500) 
end