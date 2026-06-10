function [r, p, t, l, d, e, p2, r2, p3, r3, obj] = weighted_deconv(y, delta, varargin)
%
%   This function implements the wrapper for the two-stage refinement procedure
%   described in "Artifact Aware Convex Decomposition of Electrodermal Activity for 
%   Wearable Sensors" (submitted to IEEE JBHI, June 2026)
%
%   Syntax:
%   [r, p, t, l, d, e, p2, r2, p3, r3, obj] = weighted_deconv(y, delta, alpha1, alpha2,  
%                                     alpha3, wl, psi1, psi2, psi3)
%
%   where:
%      y: observed EDA signal (we recommend to apply zscore normalization first)
%      delta: sampling interval (in seconds) of y
%      alpha1: penalization for the phasic SMNA driver (default 0.0008)
%      alpha2: penalization for the negative artifact driver (default 0.5)
%      alpha3: penalization for the positive artifact driver (default 0.5)
%      wl: duration (in seconds) of the weight window following
%          each detected artifact (default 1.3)
%      psi1: weight applied to the phasic driver penalty within artifact
%            windows in the second stage (default 0.0008)
%      psi2: weight applied to the positive artifact driver penalty within
%            artifact windows in the second stage (default 0.05)
%      psi3: weight applied to the negative artifact driver penalty within
%            artifact windows in the second stage (default 0.0005)
%
%   returns (see paper for details):
%      r: phasic component
%      p: sparse SMNA driver of phasic component
%      t: tonic component
%      l: coefficients of tonic spline
%      d: offset and slope of the linear drift term
%      e: model residuals
%      p2: sparse driver of the negative artifact component
%      r2: reconstructed negative artifact component
%      p3: sparse driver of the positive artifact component
%      r3: reconstructed positive artifact component
%      obj: value of objective function being minimized (eq 14 of paper)

% ______________________________________________________________________________
%
% File:                         weighted_deconv.m
% Last revised:                 06 June 2026 
% ______________________________________________________________________________
%
% Copyright (C) 2025-2026 Gianluca Rho, Luca Citi, Enzo Pasquale Scilingo, Alberto Greco
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
% This method is based on the established state-of-the-art cvxEDA algorithm:
% A Greco, G Valenza, A Lanata, EP Scilingo, and L Citi
% "cvxEDA: a Convex Optimization Approach to Electrodermal Activity Processing"
% IEEE Transactions on Biomedical Engineering, 2015
% DOI: 10.1109/TBME.2015.2474131
%
% This method was published in:
% G Rho, L Citi, EP Scilingo, and A Greco
% "Artifact Aware Convex Decomposition of Electrodermal Activity for Wearable Sensors"
% Submitted to IEEE Journal of Biomedical Health and Informatics, 2026
%   
%
% If you use this program in support of published research, please include a
% citation of the references above. If you use this code in a software package,
% please explicitly inform the end users of this copyright notice and ask them
% to cite the reference above in their published research.
% ______________________________________________________________________________

% parse arguments
params = {8e-4, 5e-1, 5e-1, 1.3, 8e-4, 5e-2, 5e-4};
i = ~cellfun(@isempty, varargin);
params(i) = varargin(i);
[alpha1, alpha2, alpha3, wl, psi1, psi2, psi3] = deal(params{:});


weights = [];
[r, p, t, l, d, e, p2, r2, p3, r3, obj] = ...
        artifact_cvxEDA(y, delta, alpha2, alpha3, weights, 2, 0.7, 10, alpha1);

% artifacts posthoc processing and sparsity penalties tuning
srate = 1/delta;
peaks_neg = find(p2>.25);
peaks_pos = find(p3>.25);

% if there are no artifacts then return
if isempty(peaks_neg) && isempty(peaks_pos)
    return
end

fprintf('\n There are artifacts in the EDA signal');
fprintf('\n Performing deconvolution again ... \n');
duration = wl;

weights_neg = alpha2*ones(1, length(y));
weights_pos = alpha3*ones(1, length(y));
weights_phasic = alpha1*ones(1, length(y));
for ii = 1:length(peaks_neg)
    if peaks_neg(ii)+duration*srate > length(y)
        interval = round(peaks_neg(ii):length(y));
    else
        interval = round(peaks_neg(ii):peaks_neg(ii)+duration*srate);
    end
    weights_neg(interval) = psi3;
    weights_phasic(interval) = psi1;
end
for ii = 1:length(peaks_pos)
    if peaks_pos(ii)+duration*srate > length(y)
        interval = round(peaks_pos(ii):length(y));
    else
        interval = round(peaks_pos(ii):peaks_pos(ii)+duration*srate);
    end
    weights_pos(interval) = psi2;
    weights_phasic(interval) = psi1;
end

% phasic; neg; pos
weights = [weights_phasic; weights_neg; weights_pos];
[r, p, t, l, d, e, p2, r2, p3, r3, obj] = ...
        artifact_cvxEDA(y, delta, alpha2, alpha3, weights, 2, 0.7, 10, alpha1);


end
