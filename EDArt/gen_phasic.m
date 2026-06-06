function [phasic_vec, awgn, smna_vec] = gen_phasic(N, lambda, fs)

    N = N*fs;
    snr = 33;
    meanPhasic_log = -.496;      % log microS
    sdPhasic_log = 0.200;        % log microS

    T0 = 0;
    smna_vec = zeros(N, 1);
    phasic_vec = zeros(N, 1);
    nt = 1;
    while T0 < N
        % determine time spent wout a spike (in samples)
        T = round(exprnd(lambda)) * fs;
        if T0 == 0 && T<3*fs
            continue;
        end
        if T<1*fs         % go on if the clean signal interval lasts less than Xsec
            continue;
        end
        T0 = T + T0;
        
        if T0>=N 
            break; 
        end

        % generate phasic response with Bateman ARMA model
        tau1 = 2;%unifrnd(2, 4);
        tau0 = .7;
        a1 = 1/min(tau1, tau0); % a1 > a0
        a0 = 1/max(tau1, tau0);
        delta = 1/fs;
        ar = [(a1*delta + 2) * (a0*delta + 2), 2*a1*a0*delta^2 - 8, ...
               (a1*delta - 2) * (a0*delta - 2)] / ((a1 - a0) * delta^2);
        ma = [1 2 1];

        log_sample = normrnd(meanPhasic_log, sdPhasic_log);

        smna_vec(T0, nt) = 1/delta;
        phasic = filter(ma, ar, smna_vec(:, nt));

        % rescale phasic response
        scale = 10^(log_sample);
%         phasic = (phasic-min(phasic))/(max(phasic)-min(phasic));
        phasic_vec(:, nt) = phasic;% * scale;

        nt = nt + 1;
    end

    % additive white Gaussian noise
    a = max(max(phasic_vec));
    sigma_awgn = a^2 / 10^(snr/10);
    awgn = sigma_awgn * randn(N, 1);
end