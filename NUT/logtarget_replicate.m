function [log_target] = logtarget_replicate(xobs, Beta, tau, sigmasqalpha, nbasis)

global dimen
nBeta = nbasis + 1;
Beta_1 = zeros(nBeta,(dimen + dimen*(dimen-1)/2));
Beta_2 = zeros(nBeta,dimen*(dimen-1)/2);

nsubj = size(xobs,1);
log_target = 0;
for j=1:nsubj
   log_target = log_target + logtarget_single(xobs{j}, Beta); 
end    

if dimen==2
    log_target = log_target - Beta_1(1,1)*Beta_1(1,1)/sigmasqalpha - Beta_1(2:nBeta,1)'*Beta_1(2:nBeta,1)/tau(1)-...
                              Beta_1(1,2)*Beta_1(1,2)/sigmasqalpha - Beta_1(2:nBeta,2)'*Beta_1(2:nBeta,2)/tau(2)-...
                              Beta_1(1,3)*Beta_1(1,3)/sigmasqalpha - Beta_1(2:nBeta,3)'*Beta_1(2:nBeta,3)/tau(3)-...
                              Beta_2(1:nBeta,1)'*Beta_2(1:nBeta,1)/tau(4);         
elseif dimen==3
    log_target = log_target - Beta_1(1,1)*Beta_1(1,1)/sigmasqalpha - Beta_1(2:nBeta,1)'*Beta_1(2:nBeta,1)/tau(1)-...
                              Beta_1(1,2)*Beta_1(1,2)/sigmasqalpha - Beta_1(2:nBeta,2)'*Beta_1(2:nBeta,2)/tau(2)-...
                              Beta_1(1,3)*Beta_1(1,3)/sigmasqalpha - Beta_1(2:nBeta,3)'*Beta_1(2:nBeta,3)/tau(3)-...
                              Beta_1(1,4)*Beta_1(1,4)/sigmasqalpha - Beta_1(2:nBeta,4)'*Beta_1(2:nBeta,4)/tau(4)-...
                              Beta_1(1,5)*Beta_1(1,5)/sigmasqalpha - Beta_1(2:nBeta,5)'*Beta_1(2:nBeta,5)/tau(5)-...
                              Beta_1(1,6)*Beta_1(1,6)/sigmasqalpha - Beta_1(2:nBeta,6)'*Beta_1(2:nBeta,6)/tau(6)-...
                              Beta_2(1:nBeta,1)'*Beta_2(1:nBeta,1)/tau(7)-...    
                              Beta_2(1:nBeta,2)'*Beta_2(1:nBeta,2)/tau(8)-...
                              Beta_2(1:nBeta,3)'*Beta_2(1:nBeta,3)/tau(9);                                      
end