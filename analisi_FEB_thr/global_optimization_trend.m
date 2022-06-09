%% New Optimization approach 

% global thr mean
Vth_m = pd_011.mu;

min = 0.001;
max = 20;
incr = 0.001;

sigmas = zeros(length(min:incr:max), 1);
mus = zeros(length(min:incr:max), 1);

count = 1;

for Vth_gb = min:incr:max
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
        
        Thr_no_out = Thr_best(iasic+1).mat(~isoutlier(Thr_011(iasic+1).mat(:, 7)),7);
        
        %add sample to global arrays
        %Thr_no_out = rmoutliers(Thr_no_out);
        count_opt_arr(iasic) = length(Thr_no_out);
        asics_all_channels = [asics_all_channels;Thr_no_out];
        asics_tau_means_kev(iasic) = nanmean(Thr_no_out);
    
    end
    
    Vth_ma_kev = asics_tau_means_kev;
    %Vth_ma_kev = rmoutliers(Vth_ma_kev);
    
    % global thr mean
    %Vth_m = nanmean(Vth_ma_kev)
    
    % differences
    D_kev = Vth_ma_kev - Vth_m;
    flag_counter = 0;
    changed_1 = true;
    changed_2 = true;
    
    while (changed_1 || changed_2)
    
        changed_1 = false;
        changed_2 = false;
    
        for i = 1:length(D_kev)
            if D_kev(i) > Vth_gb
                Vth_ma_kev(i) = Vth_ma_kev(i) - Vth_gb;
                corrections_ASIC(i, 1) = corrections_ASIC(i, 1) + 1;
                changed_1 = true;
            elseif D_kev(i) < -Vth_gb
                Vth_ma_kev(i) = Vth_ma_kev(i) + Vth_gb; 
                corrections_ASIC(i, 2) = corrections_ASIC(i, 2) + 1;
                changed_2 = true;
            end
        end
    
        flag_counter = flag_counter + 1;
        D_kev = Vth_ma_kev - Vth_m;
    end

    final_dist_1 = fitdist(asics_all_channels, 'Normal');
    
    
    % implement correction on channels
    
    acc_opt_temp = asics_all_channels;
    start = 1;
    
    for iasic=gbl_thr_range
        stop = start + count_opt_arr(iasic)-1;
        acc_opt_temp(start:stop) = acc_opt_temp(start:stop) - Vth_gb * corrections_ASIC(iasic, 1);
        acc_opt_temp(start:stop) = acc_opt_temp(start:stop) + Vth_gb * corrections_ASIC(iasic, 2);
    
        start = stop + 1;
    end

    final_dist_2 = fitdist(acc_opt_temp, 'Normal');
    
    mus(count, 1) = final_dist_2.mu - final_dist_1.mu;
    sigmas(count, 1) = final_dist_2.sigma - final_dist_1.sigma;

    count = count + 1;
end

figure
plot(min:incr:max, sigmas, min:incr:max, mus)
legend("\Delta\mu", "\Delta\sigma")