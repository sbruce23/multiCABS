function [Beta_out, m_out, m, xobs_tmp] = Hamilt(j, xobs, tau_temp, Beta_temp, xi_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does HMC updates for the coefficents that are different
%
%   Input:
%       1) j - one of two segments
%       2) yobs - time series observations 
%       3) tau_temp - current smoothing parameters
%       4) Beta_temp - current coefficients
%       5) xi_temp - current partitions
%   Main Outputs:
%       1) Beta_out - vector of coefficients 
%       2) m_out - updated momentum variables
%       3) m - initial momentum variables
%       4) yobs_tmp - time series observations in selected segment
%
%   Required programs:  Gradient1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    global nbasis sigmasqalpha dimen nBeta M ee

    %pick right portion of the data
    %pick right portion of the data
    if j>1
            xobs_tmp = xobs((xi_temp(j-1)+1):xi_temp(j));
    else
            xobs_tmp = xobs(1:xi_temp(j));
    end
    
    select = 1:dimen^2;
    
    ll = nBeta*length(select);
    [gr]=gradient_replicate(xobs_tmp,tau_temp, Beta_temp, sigmasqalpha, nbasis);
    Beta_old = reshape(Beta_temp(:,select),ll,1); 
    Beta_out = Beta_old;
    
    % generate momentum variable
    m = mvnrnd(zeros(ll,1),M*eye(ll))';
    % determine leap number and step
    stepsize = unifrnd(0,2*ee);
    leap = randsample(1:2*ceil(1/ee),1);
    %leap = randsample(1:(1/stepsize),1);
    
    m_out = m + 0.5*gr*stepsize;
    for i=1:leap
        Beta_out = Beta_out + stepsize*(1/M)*eye(ll)*m_out;
        Beta_temp(:,select) = reshape(Beta_out,nBeta,length(select));
        if i==leap
            [gr] = gradient_replicate(xobs_tmp, tau_temp, Beta_temp, sigmasqalpha, nbasis);
            m_out = m_out + 0.5*gr*stepsize;
        else
            [gr] = gradient_replicate(xobs_tmp, tau_temp, Beta_temp, sigmasqalpha, nbasis);
            m_out = m_out + 1*gr*stepsize;
        end    
    end    
    m_out = -m_out;
  
