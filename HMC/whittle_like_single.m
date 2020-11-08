function [log_whittle] = whittle_like_single(yobs, Beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate local Whittle likelihood
%
%   Input:
%       1) yobs - single multivariate time series
%       2) Beta - coefficient of basis functions
%   Main Outputs:
%       1) log_whitle - log local whittle likelihood
%   Require programs: lin_basis_function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dimen
Beta_1 = Beta(:,1:(dimen + dimen*(dimen-1)/2));
Beta_2 = Beta(:,(dimen + dimen*(dimen-1)/2 + 1): end);

dim = size(yobs);
n = dim(1);
nfreq = floor(n/2);
tt = (0:nfreq)/(2*nfreq);
yy = fft(yobs)/sqrt(n);
y = yy(1:(nfreq+1),:);
nf = length(y);

[xx_r, xx_i]=lin_basis_func(tt);

%theta's
theta = zeros(dimen*(dimen-1)/2,nf);
for i=1:(dimen*(dimen-1)/2)
    theta_real = xx_r * Beta_1(:,i+dimen);
    theta_imag = xx_i * Beta_2(:,i);
    theta(i,:) = theta_real + sqrt(-1)*theta_imag;
end    

%delta's
delta_sq = zeros(dimen,nf);
for i=1:dimen
    delta_sq(i,:) = exp(xx_r * Beta_1(:,i));
end    

if dimen==2
    if (mod(n,2)==1) %odd n
           log_whittle = -sum(log(delta_sq(1,2:end))' + log(delta_sq(2,2:end))' + ...
                        conj(y(2:end,1)).*y(2:end,1).*exp(-xx_r(2:end,:)*Beta_1(:,1)) + ...
                        conj(y(2:end,2) - theta(2:end).'.*y(2:end,1)).*(y(2:end,2) - theta(2:end).'.*y(2:end,1)).*exp(-xx_r(2:end,:)*Beta_1(:,2))) - ...
                        0.5*(log(delta_sq(1,1))' + log(delta_sq(2,1))' + conj(y(1,1)).*y(1,1).*exp(-xx_r(1,:)*Beta_1(:,1)) - ...
                        conj(y(1,2) - theta(1).'.*y(1,1)).*(y(1,2) - theta(1).'.*y(1,1)).*exp(-xx_r(1,:)*Beta_1(:,2)));
    else
           log_whittle = -sum(log(delta_sq(1,2:nfreq))' + log(delta_sq(2,2:nfreq))' + ...
                        conj(y(2:nfreq,1)).*y(2:nfreq,1).*exp(-xx_r(2:nfreq,:)*Beta_1(:,1)) + ...
                        conj(y(2:nfreq,2) - theta(2:nfreq).'.*y(2:nfreq,1)).*(y(2:nfreq,2) - theta(2:nfreq).'.*y(2:nfreq,1)).*exp(-xx_r(2:nfreq,:)*Beta_1(:,2))) - ...
                        0.5*(log(delta_sq(1,1)) + log(delta_sq(2,1)) + conj(y(1,1)).*y(1,1).*exp(-xx_r(1,:)*Beta_1(:,1)) + ...
                        conj(y(1,2) - theta(1).'.*y(1,1)).*(y(1,2) - theta(1).'.*y(1,1)).*exp(-xx_r(1,:)*Beta_1(:,2))) - ...
                        0.5*(log(delta_sq(1,end)) + log(delta_sq(2,end)) + conj(y(end,1)).*y((nfreq+1),1).*exp(-xx_r(end,:)*Beta_1(:,1)) + ...
                        conj(y(end,2) - theta(end).'.*y(end,1)).*(y(end,2) - theta(end).'.*y(end,1)).*exp(-xx_r(end,:)*Beta_1(:,2)));
    end
elseif dimen==3
    if (mod(n,2)==1) %odd n
        log_whittle = -sum(log(delta_sq(1,2:end))' + log(delta_sq(2,2:end))' + log(delta_sq(3,2:end))' + ...
                        conj(y(2:end,1)).*y(2:end,1).*exp(-xx_r(2:end,:)*Beta_1(:,1)) + ...
                        conj(y(2:end,2) - theta(1,2:end).'.*y(2:end,1)).*(y(2:end,2) - theta(1,2:end).'.*y(2:end,1)).*exp(-xx_r(2:end,:)*Beta_1(:,2))+...
                        conj(y(2:end,3) -(theta(2,2:end).'.*y(2:end,1)+ theta(3,2:end).'.*y(2:end,2))).*(y(2:end,3) -(theta(2,2:end).'.*y(2:end,1)+ theta(3,2:end).'.*y(2:end,2))).*...
                        exp(-xx_r(2:end,:)*Beta_1(:,3))) - ...
                        0.5*(log(delta_sq(1,1))' + log(delta_sq(2,1))' + log(delta_sq(3,1))' + ...
                        conj(y(1,1)).*y(1,1).*exp(-xx_r(1,:)*Beta_1(:,1)) + ...
                        conj(y(1,2) - theta(1,1).'.*y(1,1)).*(y(1,2) - theta(1,1).'.*y(1,1)).*exp(-xx_r(1,:)*Beta_1(:,2))+...
                        conj(y(1,3) -(theta(2,1).'.*y(1,1)+ theta(3,1).'.*y(1,2))).*(y(1,3) -(theta(2,1).'.*y(1,1)+ theta(3,1).'.*y(1,2))).*...
                        exp(-xx_r(1,:)*Beta_1(:,3)));

    else
        log_whittle = -sum(log(delta_sq(1,2:nfreq))' + log(delta_sq(2,2:nfreq))' + log(delta_sq(3,2:nfreq))' + ...
                        conj(y(2:nfreq,1)).*y(2:nfreq,1).*exp(-xx_r(2:nfreq,:)*Beta_1(:,1)) + ...
                        conj(y(2:nfreq,2) - theta(1,2:nfreq).'.*y(2:nfreq,1)).*(y(2:nfreq,2) - theta(1,2:nfreq).'.*y(2:nfreq,1)).*exp(-xx_r(2:nfreq,:)*Beta_1(:,2))+...
                        conj(y(2:nfreq,3) -(theta(2,2:nfreq).'.*y(2:nfreq,1)+ theta(3,2:nfreq).'.*y(2:nfreq,2))).*(y(2:nfreq,3) -(theta(2,2:nfreq).'.*y(2:nfreq,1)+ theta(3,2:nfreq).'.*y(2:nfreq,2))).*...
                        exp(-xx_r(2:nfreq,:)*Beta_1(:,3))) - ...
                        0.5*(log(delta_sq(1,1))' + log(delta_sq(2,1))' + log(delta_sq(3,1))' + ...
                        conj(y(1,1)).*y(1,1).*exp(-xx_r(1,:)*Beta_1(:,1)) + ...
                        conj(y(1,2) - theta(1,1).'.*y(1,1)).*(y(1,2) - theta(1,1).'.*y(1,1)).*exp(-xx_r(1,:)*Beta_1(:,2))+...
                        conj(y(1,3) -(theta(2,1).'.*y(1,1)+ theta(3,1).'.*y(1,2))).*(y(1,3) -(theta(2,1).'.*y(1,1)+ theta(3,1).'.*y(1,2))).*...
                        exp(-xx_r(1,:)*Beta_1(:,3)))-...
                        0.5*(log(delta_sq(1,end))' + log(delta_sq(2,end))' + log(delta_sq(3,end))' + ...
                        conj(y(end,1)).*y(end,1).*exp(-xx_r(end,:)*Beta_1(:,1)) + ...
                        conj(y(end,2) - theta(1,end).'.*y(end,1)).*(y(end,2) - theta(1,end).'.*y(end,1)).*exp(-xx_r(end,:)*Beta_1(:,2))+...
                        conj(y(end,3) -(theta(2,end).'.*y(end,1)+ theta(3,end).'.*y(end,2))).*(y(end,3) -(theta(2,end).'.*y(end,1)+ theta(3,end).'.*y(end,2))).*...
                        exp(-xx_r(end,:)*Beta_1(:,3)));
                
    end
end  
%log_whittle = log_whittle/n;