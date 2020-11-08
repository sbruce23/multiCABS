function [Beta_out,xobs_tmp] = NUT_sampler(j, xobs, tau_temp, Beta_temp, xi_temp, uu_cut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does HMC updates for the coefficents that are different
%
%   Input:
%       1) j - one of two segments
%       2) yobs - a cell contains N replicated P dimensional time series
%       3) tau_temp - current smoothing parameters
%       4) Beta_temp - current coefficients
%       5) xi_temp - current partitions
%       6?uu_cut - group partitions
%   Main Outputs:
%       1) Beta_out - vector of coefficients 
%       2) m_out - updated momentum variables
%       3) m - initial momentum variables
%       4) yobs_tmp - time series observations in selected group
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    global nbasis sigmasqalpha dimen nBeta

    %pick right portion of the data
    if j>1
            xobs_tmp = xobs((uu_cut(xi_temp(j-1))+1):uu_cut(xi_temp(j)));
    else
            xobs_tmp = xobs(1:uu_cut(xi_temp(j)));
    end
    
    
    f = @(Beta_temp) deal(logtarget_replicate(xobs_tmp, Beta_temp, tau_temp, sigmasqalpha, nbasis),...
                gradient_replicate(xobs_tmp, Beta_temp, tau_temp, sigmasqalpha, nbasis));
    Beta_old = reshape(Beta_temp,nBeta*dimen^2,1);        
    %[Beta_tmp, stepsize, ~, ~] = dualAveraging(f, Beta_old, 0.8, 10);
    [Beta_out, ~, ~, ~] = NUTS(f, 0.1, Beta_old);