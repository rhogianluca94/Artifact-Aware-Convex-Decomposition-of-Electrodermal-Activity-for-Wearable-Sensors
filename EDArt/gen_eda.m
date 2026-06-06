function [y, t, r, p, art, artTimes] = gen_eda(lambda_scr, lambda_art, L, fs)

    % Inputs 
    % lambda_scr: expected time occurring between two consecutive SCRs
    % lambda_art: expected time occurring between two artifacts
    % L: length of generated EDA in seconds
    % fs: sampling frequency
    %
    % Outputs
    % y: simulated EDA
    % t: tonic component
    % r: phasic component
    % p: SMNA neural bursts (unitary amplitude)
    % art: hand-movement artifacts

    % phasic component
    [phasic_vec, noise, smna_vec] = gen_phasic(L, lambda_scr, fs);
    phasic = sum(phasic_vec, 2);
    p = sum(smna_vec, 2);

%     fc = 0.08;   % Cutoff frequency (Hz)
%     order = 1;   % Filter order 
%     Rp = 0.5;    % Passband ripple (dB)
%     Rs = 80;     % Stopband attenuation (dB)
%     [b, a, k] = ellip(order, Rp, Rs, 2*fc/fs, 'low');
%     sos = zp2sos(b, a, k); 
%     r = sosfilt(sos, phasic);
    r = phasic;

    % tonic component
    t = gen_tonic(L, fs);

    % artifacts component
    [artifact_vec, artTimes, artTypes] = gen_prob_states(lambda_art, L, fs);
    art = sum(artifact_vec, 2);

    y = t + r + art + noise;

end