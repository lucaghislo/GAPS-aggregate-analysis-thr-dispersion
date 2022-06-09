%% Cubic gain analysis for each channel and ASIC (on FEB)

% clear; clc;

main_folder = 'aggregated_gain/';

if (~exist(main_folder,'dir'))
    mkdir(main_folder);
end

% % conversion params
% 
% conv = 31.25 %[V / DAC_inj code]
% kn = 1.6e-19 / 3.6e-3 %[C / keV]
% gamma = (conv * C_inj * 10**(-6)) / kn %[kev / DAC_inj code];
as_n_av = 0;    
N = 538;
cubic_gain = zeros(32,N);
class_A_idxs = nan(8,N);
class_B_idxs = nan(8,N);
for j=0:7
    for i=0:N-1
        number = sprintf('%03i',i);
        pathtofile = ['FEB_Test/MODULE_',number,'/1/analysis_matlab/TransferFunction/Low_Energy_Gain.dat'];
        if (exist(pathtofile))
            data = importdata(pathtofile);
            cubic_gain(:,i+1) = data.data(j+1:8:end,15);
        else
            % disp(["FEB" number "data not available in the location..."]);
            as_n_av = as_n_av + 1;
        end
    end
    
    cubic_gain(isnan(cubic_gain) | cubic_gain <= 0) = nan;

    % remove invalid FEBs
    plot_cubic_gain = cubic_gain(:,~isnan(cubic_gain(1,:)));
    asic_n = 0:1:N-1;
    plot_asic_n = asic_n(~isnan(cubic_gain(1,:)));

    % remove outliers for plotting
    % plot_cubic_gain(plot_cubic_gain<0.304)=NaN;
    % plot_cubic_gain(plot_cubic_gain>1.3)=NaN;

    % gain dispersion
    % plot_cubic_gain = abs(plot_cubic_gain - mean(plot_cubic_gain));

    figure;
    set(gcf,'Visible','off');
    set(gcf,'units','normalized','innerposition',[0 0 .9 .9]);
    heatmap(plot_cubic_gain,'XDisplayLabels',plot_asic_n,'YDisplayLabels',0:1:31);
    colormap('jet');
    xlabel('#FEB');
    ylabel('#Channel');
    caxis([0.85,1.15])
    set(gca,'FontSize',8);
    title(sprintf('Cubic gain (in x=0) for each channel in each FEB [ADCu / DAC\\_inju]. Peaking time = tau\\_%i',j));
    file_out = sprintf([main_folder, 'Cubic_gain_tau%i.pdf'],j);
    exportgraphics(gcf,file_out,'ContentType','vector');
    close(gcf);
    
    file_out = sprintf([main_folder, 'Cubic_gain_tau%i.dat'],j);
    
    fileID = fopen(file_out,'w');
    fprintf(fileID,'#---------------------------------------------------------\n');
    fprintf(fileID,'# Cubic gain (in x=0) for each channel in each FEB [ADCu / DAC_inju]. Peaking time = %i\n',j);
    fprintf(fileID,'# Rows: #Channels\n');
    fprintf(fileID,'# Columns: #FEB\n');
    fprintf(fileID,'#---------------------------------------------------------\n');
    s = repmat('%2.3f\t',1,N);
    s(end+1:end+2) = '\n';
    fprintf(fileID,s,cubic_gain');
    fclose(fileID);
    
    temp = isnan(cubic_gain);
    temp = sum(temp);
    temp_A = temp == 0;
    temp_B = temp >= 1 & temp <= 1;
    class_A_idxs(j+1,:) = temp_A;
    class_B_idxs(j+1,:) = temp_B;
end


%% ASIC Classification
classification = nan(8,N);
classification(class_A_idxs == 1) = 3;
classification(class_B_idxs == 1) = 2;
classification(isnan(classification)) = 1;

temp = sum(class_A_idxs,1);
ASIC_class_A = temp == 8;
temp = class_A_idxs | class_B_idxs;
temp = sum(temp,1);
ASIC_class_B = temp == 8;
ASIC_class_B(ASIC_class_A) = 0;

figure('units','normalized','outerposition',[0 0 .9 .9], 'visible', 'off')
%subaxis(1,1,1,'Padding', 0.04, 'Margin', 0.01);
heatmap(classification,'XDisplayLabels',0:1:N-1,'YDisplayLabels',0:1:7);
xlabel('#ASIC');
ylabel('#Channel');
title(sprintf('FEB Classification. Class A = %d, Class B = %d',sum(ASIC_class_A),sum(ASIC_class_B)));
colormap([0.95,0,0; ...
          0.95,0.95,0; ...
          0,0.95,0])
exportgraphics(gcf,sprintf([main_folder, 'Classification.pdf']));
close(gcf);

tempA = find(ASIC_class_A == 1) - 1; % index starts from 1
tempB = find(ASIC_class_B == 1) - 1;
display(['Class A FEBs: ', sprintf(repmat('%i\t',1,length(tempA)), tempA)])
display(['Class A FEBs total: ', sprintf('%d', sum(ASIC_class_A))])
display(['Class B FEBs: ', sprintf(repmat('%i\t',1,length(tempB)), tempB)])
display(['Class B FEBs total: ', sprintf('%d', sum(ASIC_class_B))])
display(['FEBs not available: ', sprintf('%d', as_n_av/8)])

file_out = [main_folder, 'FEB_Classification_A.dat'];

fileID = fopen(file_out,'w');
fprintf(fileID,'#---------------------------------------------------------\n');
fprintf(fileID,'# A-class FEBs.\n');
fprintf(fileID,'# Columns: #FEB\n');
fprintf(fileID,'#---------------------------------------------------------\n');
s = repmat('%i\t',1,length(tempA));
fprintf(fileID,s,tempA);
fclose(fileID);

file_out = [main_folder, 'FEB_Classification_B.dat'];

fileID = fopen(file_out,'w');
fprintf(fileID,'#---------------------------------------------------------\n');
fprintf(fileID,'# B-class FEBs.\n');
fprintf(fileID,'# Columns: #FEB\n');
fprintf(fileID,'#---------------------------------------------------------\n');
s = repmat('%i\t',1,length(tempB));
fprintf(fileID,s,tempB);
fclose(fileID);
