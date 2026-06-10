function [r, p, t, l, d, e, p2, r2, p3, r3, obj] = artifact_cvxEDA(y, delta, alpha2, alpha3, weights, varargin)
%
%   This function implements the core convex optimization algorithm described
%   in "Artifact Aware Convex Decomposition of Electrodermal Activity for 
%   Wearable Sensors" (submitted to IEEE JBHI, June 2026)
%
%   Syntax:
%   [r, p, t, l, d, e, obj] = artifact_cvxEDA(y, delta, alpha2, alpha3, weights,  
%                                    tau0, tau1, delta_knot, alpha, gamma)
%
%   where:
%      y: observed EDA signal (we recommend normalizing it: y = zscore(y))
%      delta: sampling interval (in seconds) of y
%      tau0: slow time constant of the Bateman function (default 2.0)
%      tau1: fast time constant of the Bateman function (default 0.7)
%      delta_knot: time between knots of the tonic spline function (default 10)
%      alpha: penalization for the sparse SMNA driver (default 0.0008)
%      gamma: penalization for the tonic spline coefficients (default 0.01)
%
%   returns (see paper for details):
%      r: phasic component
%      p: sparse SMNA driver of phasic component
%      t: tonic component
%      l: coefficients of tonic spline
%      d: offset and slope of the linear drift term
%      e: model residuals
%      obj: value of objective function being minimized (eq 14 of paper)

% ______________________________________________________________________________
%
% File:                         artifact_cvxEDA.m
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
% This method was first published in:
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
params = {2, 0.7, 10, 8e-4, 1e-2};
i = ~cellfun(@isempty, varargin);
params(i) = varargin(i);
[tau0, tau1, delta_knot, alpha, gamma] = deal(params{:});

n = length(y);
y = y(:);

% bateman ARMA model
a1 = 1/min(tau1, tau0); % a1 > a0
a0 = 1/max(tau1, tau0);
ar = [(a1*delta + 2) * (a0*delta + 2), 2*a1*a0*delta^2 - 8, ...
       (a1*delta - 2) * (a0*delta - 2)] / ((a1 - a0) * delta^2);
ma = [1 2 1];

% matrices for ARMA model
i = 3:n;
A = sparse([i i i], [i i-1 i-2], repmat(ar, n-2, 1), n, n);
M = sparse([i i i], [i i-1 i-2], repmat(ma, n-2, 1), n, n);

b = (3.3680) * 2.4481;
a = (3.8202) * 2.3378;
ampl = 2.3542 * .1919;
c1 = exp(-a*delta)*delta;
c2 = exp(-b*delta)*delta;
ar_art = [1, ...
      -2*exp(-b*delta)-2*exp(-a*delta), ...
      4*exp(-(a+b)*delta) + exp(-2*b*delta)+exp(-2*a*delta),...
      -2*exp(-(a+2*b)*delta)-2*exp(-(2*a+b)*delta), ...
      exp(-2*(a+b)*delta)];
ma_art = [1 ...
      c1+c2 ...
      -2*c1*exp(-b*delta)-2*c2*exp(-a*delta) ...
      exp(-2*b*delta) + exp(-2*a*delta)] ./(2000) * ampl * 35;%./ 3636.11* ampl;

% positive artifact
i = 5:n;
A_art = sparse([i i i i i], [i i-1 i-2 i-3 i-4], repmat(ar_art, n-4, 1), n, n);
i = 4:n;
M_art = sparse([i i i i], [i i-1 i-2 i-3], repmat(-ma_art, n-3, 1), n, n);

% negative artifact
M_art_p = - M_art;
A_art_p = A_art;

% spline
delta_knot_s = round(delta_knot / delta);
spl = [1:delta_knot_s delta_knot_s-1:-1:1]'; % order 1
spl = conv(spl, spl, 'full');
spl = spl / max(spl);
% matrix of spline regressors
i = bsxfun(@plus, (0:length(spl)-1)'-floor(length(spl)/2), 1:delta_knot_s:n);
nB = size(i, 2);
j = repmat(1:nB, length(spl), 1);
p = repmat(spl(:), 1, nB);
valid = i >= 1 & i <= n;
B = sparse(i(valid), j(valid), p(valid));

% trend
C = [ones(n,1) (1:n)'/n];
nC = size(C, 2);

% Solve the problem:
% .5*(M*q + B*l + C*d + M_art*q_art + M_art_p*q_art_p - y)^2 + 
% alpha*sum(A,1)*q + alpha2*sum(A_art,1)*q_art  + alpha3*sum(A_art_p,1)*q_art_p + .5*gamma*l'*l
% s.t. A*q >= 0
%      A_art*q_art >= 0 

% two artifacts
H = [M'*M, M'*C, M'*B M'*M_art, M'*M_art_p; ...
    C'*M, C'*C, C'*B, C'*M_art, C'*M_art_p; ...
    B'*M, B'*C, B'*B+gamma*speye(nB), B'*M_art, B'*M_art_p; ...
    M_art'*M, M_art'*C, M_art'*B, M_art'*M_art, M_art'*M_art_p; ...
    M_art_p'*M, M_art_p'*C, M_art_p'*B, M_art_p'*M_art, M_art_p'*M_art_p];

if isempty(weights)
    % scalar weighting of sparsity penalty parameters
    term_artNeg = alpha3 * sum(A_art,1)';
    term_artPos = alpha2 * sum(A_art_p,1)';
    term_phasic = alpha * sum(A, 1)';
else
    % point-wise weighting of sparsity penalty parameters
    term_phasic = sum(A, 1)' .* weights(1,:)';
    term_artNeg = sum(A_art,1)' .* weights(2,:)';
    term_artPos = sum(A_art_p,1)' .* weights(3,:)';
end

f = [term_phasic-M'*(y); -(C'*(y)); -(B'*(y)); term_artNeg-M_art'*y; term_artPos-M_art_p'*y];

H = H + 1e-8*speye(size(H));

N = 3*n+nB+nC; % coeffs
Aineq = sparse([-A zeros(n, N-n); ...
                zeros(n, N-2*n) -A_art zeros(n, n); ...
                zeros(n, N-n) -A_art_p]);
bineq = sparse(zeros(3*n, 1));


options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', ...
    'Display', 'final', 'ConstraintTolerance', 1e-12, 'OptimalityTolerance', 1e-12, 'MaxIterations', 500);
[z, obj] = quadprog(H, f, Aineq, bineq, ...
    [], [], [], [], [], options);
obj = obj + .5 * (y' * y);

q = z(1:n); 
d = z(n+1:n+nC);
l = z(n+nC+1:n+nB+nC);
q_art = z(n+nB+nC+1:end-n);
q_art_p = z(end-n+1:end);

p = A * q;
r = M * q;
p2 = A_art * q_art;
r2 = M_art * q_art;
p3 = A_art_p * q_art_p;
r3 = M_art_p * q_art_p;
t = B*l + C*d;
e = y - r - t - r2 - r3;

end
