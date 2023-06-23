function [val,idx] = of_NRMSE(obs,sim,idx)
% of_NRMSE Calculates the Normalized Root Mean Squared Error of simulated 
% streamflow, with normalization by the range of observations. Ignores time
% steps with negative flow values.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% idx       - optional vector of indices to use for calculation, can be
%               logical vector [nx1] or numeric vector [mx1], with m <= n
%
% Out:
% val       - objective function value          [1x1]
% idx       - vector of indices used for calculation


%% Check inputs and select timesteps
if nargin < 2
    error('Not enugh input arguments')    
end

if nargin < 3; idx = []; end
[sim, obs, idx] = check_and_select(sim, obs, idx);                                             

%% Calculate metric
norm = max(obs)-min(obs);
val = sqrt(mean((obs-sim).^2)) / norm;
end

