function[Beta_mean,Beta_var,xobs_tmp] = postBeta(j, xobs, tau_temp, xi_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the mean and variance of normal approximation for coefficient
% of basis fucntions that are different across segments
%
%   Input:
%       1) j - one of two segments
%       2) yobs - time series observations 
%       3) tau_temp - current smoothing parameters
%       4) Beta_temp - current coefficients
%       5) xi__temp - current partitions
%   Main Outputs:
%       1) Beta_mean - mean of approximated distribution
%       2) Beta_var - variance of approximated distribution
%       3) yobs_tmp - time series observations in selected segment
%
%   Required programs:  Beta_derive1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    global nbasis nBeta sigmasqalpha dimen options
   
    
    %pick right portion of the data
    if j>1
            xobs_tmp = xobs((xi_temp(j-1)+1):xi_temp(j));
    else
            xobs_tmp = xobs(1:xi_temp(j));
    end
    
    %provide initial values
    x = zeros(nBeta*dimen^2,1);
    
    %optimization process
    [Beta_mean,~,~,~,~,Beta_inv_var] = fminunc(@Beta_derive_replicate, x, options, ...
                        xobs_tmp, tau_temp, sigmasqalpha, nbasis); 
    %Beta_var = Beta_inv_var\eye(size(Beta_inv_var));
    Beta_var = matpower(Beta_inv_var,-1);
    