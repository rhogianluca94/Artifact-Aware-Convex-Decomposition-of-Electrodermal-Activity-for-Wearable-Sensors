function tonic = gen_tonic(N, fs)
    % Generate tonic component filtering white noise with an IIR filter

    fc = 0.01;   % Cutoff frequency (Hz)
    order = 6;   % Filter order 
    Rp = 0.5;    % Passband ripple (dB)
    Rs = 80;     % Stopband attenuation (dB)
    [z, p, k] = ellip(order, Rp, Rs, 2*fc/fs, 'low');
    sos = zp2sos(z, p, k); % Convert to Second-Order Sections
    
    meanTonic_log = 2.174;      % log kOhm
    sdTonic_log = .1;%0.205;        % log kOhm
    tonic_max_range = 40;       % microS
    tonic_min_range = 1;        % microS
    log_sample = normrnd(meanTonic_log, sdTonic_log);
    offset = (1 / 10^log_sample) * 1e3;  % Convert to microS
    scale = unifrnd(tonic_min_range, min(offset-tonic_min_range, tonic_max_range));
    wn = cumsum(randn(1, N*fs));
    tonic = sosfilt(sos, wn)';
    tonic = (tonic-min(tonic))/(max(tonic)-min(tonic));
    tonic = tonic*scale + offset;
end