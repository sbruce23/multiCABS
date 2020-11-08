function [spect]=MAspectriv(beta,sigma,freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take trivariate MA coefficient matrix, output the spectral
% beta: 3x(3xlag) coefficient matrix
% sigma: 3x3 covariance matrix
% freq: freqency where spectral calcuated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


phispect = zeros(3,3,length(freq));
spect = zeros(3,3,length(freq));

for k=1:length(freq)
    phispect(:,:,k) = eye(3);
    for j=1:2
        if j==1
            bigmat = beta(:,1:3).* exp(-2*pi*sqrt(-1)*freq(k));
        else
            bigmat = beta(:,(3*j-2):3*j) .* exp(-2*j*pi*sqrt(-1)*freq(k));
        end
        phispect(:,:,k) = phispect(:,:,k) + bigmat;
    end
    spect(:,:,k) = (phispect(:,:,k))*sigma*(phispect(:,:,k)');
end    
