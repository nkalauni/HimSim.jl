function [out] = rainfall_1(In,T,p1,varargin)
%rainfall_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Rainfall based on temperature threshold
% Constraints:  -
% @(Inputs):    p1   - temperature threshold above which rainfall occurs [oC]
%               T    - current temperature [oC]
%               In   - incoming precipitation flux [mm/d]
%               varargin(1) - smoothing variable r (default 0.01)

if size(varargin,2) == 0
    out = In.*(1-smoothThreshold_temperature_logistic(T,p1));
elseif size(varargin,2) == 1
    out = In.*(1-smoothThreshold_temperature_logistic(T,p1,varargin(1)));    
end

end

