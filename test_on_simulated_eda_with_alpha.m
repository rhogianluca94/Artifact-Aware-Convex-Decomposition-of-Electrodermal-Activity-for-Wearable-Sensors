% This exemplary code refers to Section III-A of the paper "Artifact Aware 
% Convex Decomposition of Electrodermal Activity for Wearable Sensors" 
% (submitted to IEEE JBHI, June 2026)
%
% ______________________________________________________________________________
%
% File:                         test_on_simulated_EDA_with_alpha.m
% Last revised:                 06 June 2026 
% ______________________________________________________________________________
%
% Copyright (C) 2025-2026 Gianluca Rho, Luca Citi, Alberto Greco
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You may contact the author by e-mail (gianluca.rho@ing.unipi.it).
% ______________________________________________________________________________
%
%
% This method used in this exemplary script is based on the established 
% state-of-the-art cvxEDA algorithm:
% A Greco, G Valenza, A Lanata, EP Scilingo, and L Citi
% "cvxEDA: a Convex Optimization Approach to Electrodermal Activity Processing"
% IEEE Transactions on Biomedical Engineering, 2015
% DOI: 10.1109/TBME.2015.2474131
%
% The method tested in this exemplary script was published in:
% G Rho, L Citi, EP Scilingo, and A Greco
% "Artifact Aware Convex Decomposition of Electrodermal Activity for Wearable Sensors"
% Submitted to IEEE Journal of Biomedical Health and Informatics, 2026
%
% The method used to generate artifical artifact activity (EDArt) was 
% published in:
% G Rho, N Carbonaro, M Laurino, A Tognetti, and A Greco
% "EDArt: A Framework for the Simulation of Hand-Movement Artifact-Corrupted 
% Electrodermal Activity Signal"
% 2025 IEEE International Conference on Metrology for eXtended Reality, 
% Artificial Intelligence and Neural Engineering (MetroXRAINE), 
% pp. 1325-1330, 2025
% DOI: 10.1109/MetroXRAINE66377.2025.11340365
%   
%
% If you use this program in support of published research, please include a
% citation of the references above. If you use this code in a software package,
% please explicitly inform the end users of this copyright notice and ask them
% to cite the reference above in their published research.
% ______________________________________________________________________________
warning('off', 'signal:findpeaks:largeMinPeakHeight');

clc
close all
clear
basePath = pwd;
psave = fullfile(basePath, 'performance_results/');
addpath([basePath filesep 'EDArt'])

if ~exist(psave)
    mkdir(psave)
end

plot_results = 1;

srate = 50;
delta_knot = 10;
time = 0:1/srate:180-1/srate;
nRep = 30;
Nsurr = 1;
sigmas = [8.2e-3, 7.35e-2];
snrLevel = {'high', 'low'};

tau1 = 2;
tau0 = .7;
a1 = 1/min(tau1, tau0);
a0 = 1/max(tau1, tau0);
delta = 1/srate;
ar = [(a1*delta + 2) * (a0*delta + 2), 2*a1*a0*delta^2 - 8, ...
       (a1*delta - 2) * (a0*delta - 2)] / ((a1 - a0) * delta^2);
ma = [1 2 1];

s = struct();

% === setup accumulators for each SNR level ===
for jj = 1:length(sigmas)
    PPV_std_all{jj} = [];
    PPV_corr_all{jj} = [];
    SENS_std_all{jj} = [];
    SENS_corr_all{jj} = [];
    PPV_art_all{jj} = [];
    SENS_art_all{jj} = [];
    msePhasicCorrected_all{jj} = [];
    msePhasicStandard_all{jj} = [];
    mseTonicCorrected_all{jj} = [];
    mseTonicStandard_all{jj} = [];
end

% ================= Testing loop =================
for kk = 1:Nsurr

    % ====== Create Phasic component ======
    smnas = zeros(length(time), nRep);
    pos = [];
    while length(pos) < nRep
        posI = randi([5*srate, length(time)-5*srate]);
        if isempty(pos) || (~any(abs(pos-posI))) < 1*srate
            smnas(posI, length(pos)+1) = 1/delta;
            pos = [pos posI];
        end
    end
    smna = sum(smnas, 2);
    phasic = filter(ma, ar, smna);

    % ====== Create Tonic component ======
    sin_min = 45;
    sin_max = 90;
    tonic = 2 + rand()*linspace(0,1,length(time))' + ...
        sin(2*pi*(rand()+time/(sin_min + rand()*(sin_max - sin_min))))';

    % ====== Create Artifacts ======
    lambda_scr = 5;
    lambda_art = 10;
    L = length(tonic);
    [y, ytonic, yphasic, ysmna, yart, artp] = ...
        gen_eda(lambda_scr, lambda_art, L/srate, srate);

    for jj = 1:length(sigmas)

        % ====== noise ======
        sigma = sigmas(jj);
        noise = sigma * randn(size(y));
        y = yphasic + tonic + noise;
        yartifacted = zscore(y + yart);

        % ====== Weighted deconv ======
        alpha = 8e-4;
        alpha2 = 5e-1;
        alpha3 = 5e-1;

        [r, p, t, l, d, e, p2, r2, p3, r3, obj] = ...
            weighted_deconv(yartifacted, delta, alpha, alpha2, alpha3, []);

        % ====== cvxEDA ======
        [r_orig, p_orig, t_orig, l, d, e] = ...
            cvxEDA(zscore(y), delta, 2, 0.7, 10, alpha);

        % ====== Compute SMNA PPV and Sens ======
        [PPV_std, SENS_std] = computePPVfromSMNA(ysmna, p_orig, time);
        [PPV_corr, SENS_corr] = computePPVfromSMNA(ysmna, p, time);

        PPV_std_all{jj}  = [PPV_std_all{jj}, PPV_std];
        PPV_corr_all{jj} = [PPV_corr_all{jj}, PPV_corr];
        SENS_std_all{jj} = [SENS_std_all{jj}, SENS_std];
        SENS_corr_all{jj}= [SENS_corr_all{jj}, SENS_corr];

        % ====== Compute Artifact Drivers metrics ======
        smnaArt = zeros(size(p2));
        smnaArt(artp) = 50;
        smnaArt_est = (p2 + p3);
        [PPV_art, SENS_art] = ...
            computePPVfromSMNA(smnaArt, smnaArt_est*1/delta, time);

        PPV_art_all{jj}  = [PPV_art_all{jj}, PPV_art];
        SENS_art_all{jj} = [SENS_art_all{jj}, SENS_art];

        % ====== Mean Squared Error ======
        winLenPhasic = 5*srate;
        winLenTonic = 20*srate;

        msePhasicCorrected = [];
        msePhasicStandard = [];
        mseTonicCorrected = [];
        mseTonicStandard = [];

        nnT = 1;
        nnR = 1;
        wStartTonic = 1;

        for ii = 1:winLenPhasic:L
            wstart = ii;
            wstop = wstart + winLenPhasic - 1;

            phasicTrue = yphasic(wstart:wstop);
            msePhasicCorrected(nnR) = mean((r(wstart:wstop)-phasicTrue).^2);
            msePhasicStandard(nnR)  = mean((r_orig(wstart:wstop)-phasicTrue).^2);

            if mod(wstop-wStartTonic, winLenTonic-1) == 0
                tonicTrue = tonic(wStartTonic:wstop);
                mseTonicCorrected(nnT) = mean((t(wStartTonic:wstop)-tonicTrue).^2);
                mseTonicStandard(nnT)  = mean((t_orig(wStartTonic:wstop)-tonicTrue).^2);
                wStartTonic = wstart;
                nnT = nnT + 1;
            end
            nnR = nnR + 1;
        end

        msePhasicStandard_all{jj}  = [msePhasicStandard_all{jj}, msePhasicStandard];
        msePhasicCorrected_all{jj} = [msePhasicCorrected_all{jj}, msePhasicCorrected];
        mseTonicStandard_all{jj}   = [mseTonicStandard_all{jj}, mseTonicStandard];
        mseTonicCorrected_all{jj}  = [mseTonicCorrected_all{jj}, mseTonicCorrected];

        % ====== Plots ======
        if plot_results
            figure, 
            ax1 = subplot(5,1,1); hold on,
            plot(time, zscore(y), 'LineWidth', 2), 
            plot(time, yartifacted, 'LineWidth', 2), grid on,
        
            ax2 = subplot(5,1,2); hold on
            plot(time, r, 'LineWidth', 2, DisplayName='Deconv corrected'), grid on, title('Phasic')
            plot(time, r_orig, 'LineWidth', 2, DisplayName='Deconv orig'), legend
            plot(time, yphasic, 'LineWidth', 2, DisplayName='True'), legend
        
            ax3 = subplot(5,1,3); hold on
            plot(time, t, 'LineWidth', 2, DisplayName='Deconv'), grid on, title('Tonic')
            plot(time, t_orig, 'LineWidth', 2, DisplayName='Deconv orig'), legend
        
            ax4 = subplot(5,1,4); hold on
            plot(time, p*delta, 'LineWidth', 2, DisplayName='Deconv corrected'), grid on, title('SMNA')
            plot(time, ysmna*delta, 'LineWidth', 2, DisplayName='True'), legend
            plot(time, p_orig*delta, 'LineWidth', 2, DisplayName='Deconv orig'), legend
        
            ax5 = subplot(5,1,5); hold on
            plot(time, p2, 'LineWidth', 2, DisplayName='Deconv p2'), grid on, title('Artifact')
            plot(time, p3, 'LineWidth', 2, DisplayName='Deconv p3'), legend
        
            linkaxes([ax1 ax2 ax3 ax4 ax5], 'x')
        end

    end
end

% ====== Save ======
for jj = 1:length(sigmas)
    s.phasic.mseStandard.(snrLevel{jj})  = msePhasicStandard_all{jj};
    s.phasic.mseCorrected.(snrLevel{jj}) = msePhasicCorrected_all{jj};
    s.tonic.mseStandard.(snrLevel{jj})   = mseTonicStandard_all{jj};
    s.tonic.mseCorrected.(snrLevel{jj})  = mseTonicCorrected_all{jj};
    s.PPV_corr.(snrLevel{jj}) = PPV_corr_all{jj};
    s.PPV_std.(snrLevel{jj})  = PPV_std_all{jj};
    s.SENS_corr.(snrLevel{jj})= SENS_corr_all{jj};
    s.SENS_std.(snrLevel{jj}) = SENS_std_all{jj};
end
save([psave filesep 'performance_on_simulated_data.mat'], 's');

%--------------------------------------------------------------------------
function [PPV, SENS, TP, FP, FN] = computePPVfromSMNA(smnaTrue, smnaEst, time)

    [pks, locsTrue] = findpeaks(smnaTrue/50, 'MinPeakHeight', 0.4);
    [pks, locsEst]  = findpeaks(smnaEst/50, 'MinPeakHeight', 0.4);
    smnaEstTimes  = time(locsEst); 
    smnaTrueTimes = time(locsTrue);

    TP = 0;
    FP = 0;
    FN = 0;
    
    deltaT = 0.150;
    
    matchedTrue = false(size(smnaTrueTimes));

    for kk = 1:length(smnaEstTimes)
        timeDiffs = abs(smnaTrueTimes - smnaEstTimes(kk));
        [minDiff, minIdx] = min(timeDiffs);
        
        if minDiff <= deltaT && ~matchedTrue(minIdx)
            TP = TP + 1;
            matchedTrue(minIdx) = true;
        else
            FP = FP + 1;
        end
    end

    FN = sum(~matchedTrue);

    % PPV
    if TP + FP > 0
        PPV = TP / (TP + FP) * 100;
    else
        PPV = 0;
    end

    % Sens
    if TP + FN > 0
        SENS = TP / (TP + FN) * 100;
    else
        SENS = 0;
    end

end
