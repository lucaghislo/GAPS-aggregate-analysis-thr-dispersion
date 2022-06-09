% get mean threshold
clear mean_arr
for j=0:7
    for iasic=0:483
        mean_arr(iasic+1,j+1) = nanmean(Thr_best(iasic+1).mat(~isoutlier(Thr_best(iasic+1).mat(:, j+1)), j+1));   
    end
end
mean_arr(133,:) = nan(1,8);

%get asic w/ channel w/ gain < 0.8
clear idxs
cubic_gain_temp = cubic_gain;

%pick A-class
matA = readmatrix('ASIC_Classification_A.dat');
matA = matA(4,:);

i=1;
for j=7:7
    for iasic=matA
        cubic_gain_temp(iasic+1).mat(isnan(cubic_gain_temp(iasic+1).mat) | cubic_gain_temp(iasic+1).mat <= 0.5) = nan;
        if (sum(isnan(cubic_gain_temp(iasic+1).mat(:))) > 0)
            idxs(i) = iasic;
            i = i+1;
        end
    end
end
unique(idxs)
length(unique(idxs))