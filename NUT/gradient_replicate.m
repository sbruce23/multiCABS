function [gr] = gradient_replicate(xobs, Beta_temp, tau_temp, sigmasqalpha, nbasis)

global dimen

nsubj = size(xobs,1);
gr = 0;
for j=1:nsubj
   [gr_tmp] = gradient_single(xobs{j}, Beta_temp, tau_temp, sigmasqalpha); 
   gr = gr + gr_tmp;
end    

%initilize Beta_1 and Beta_2: Beta_2 is for imaginary components
nBeta = nbasis + 1;
Beta_temp = reshape(Beta_temp,nBeta,dimen^2);
Beta_1 = zeros(nBeta,(dimen + dimen*(dimen-1)/2));
Beta_2 = zeros(nBeta,dimen*(dimen-1)/2);
Beta_1(:,:) = Beta_temp(:,1:(dimen + dimen*(dimen-1)/2));
Beta_2(:,:) = Beta_temp(1:nBeta,(dimen + dimen*(dimen-1)/2 + 1): end);

if dimen==2
      
    gr1 = zeros(nBeta,1); gr2 = zeros(nBeta,1); gr3 = zeros(nBeta,1); gr4 = zeros(nBeta,1);
    gr1(1) = Beta_1(1,1)/sigmasqalpha; gr1(2:nBeta,1) = Beta_1(2:nBeta,1)/tau_temp(1);
    gr2(1) = Beta_1(1,2)/sigmasqalpha; gr2(2:nBeta,1) = Beta_1(2:nBeta,2)/tau_temp(2);
    gr3(1) = Beta_1(1,3)/sigmasqalpha; gr3(2:nBeta,1) = Beta_1(2:nBeta,3)/tau_temp(3);
    gr4(1:nBeta,1) = Beta_2(1:nBeta,1)/tau_temp(4);   
    grs = [gr1;gr2;gr3;gr4];
    gr = gr  + (nsubj-1)*grs;
    
else
    
    gr1 = zeros(nBeta,1); gr2 = zeros(nBeta,1); gr3 = zeros(nBeta,1); gr4 = zeros(nBeta,1);
    gr5 = zeros(nBeta,1); gr6 = zeros(nBeta,1); gr7 = zeros(nBeta,1); gr8 = zeros(nBeta,1);
    gr9 = zeros(nBeta,1);

    gr1(1) = Beta_1(1,1)/sigmasqalpha; gr1(2:nBeta) = Beta_1(2:nBeta,1)/tau_temp(1);
    gr2(1) = Beta_1(1,2)/sigmasqalpha; gr2(2:nBeta) = Beta_1(2:nBeta,2)/tau_temp(2);
    gr3(1) = Beta_1(1,3)/sigmasqalpha; gr3(2:nBeta) = Beta_1(2:nBeta,3)/tau_temp(3);
    gr4(1) = Beta_1(1,4)/sigmasqalpha; gr4(2:nBeta) = Beta_1(2:nBeta,4)/tau_temp(4);
    gr5(1) = Beta_1(1,5)/sigmasqalpha; gr5(2:nBeta) = Beta_1(2:nBeta,5)/tau_temp(5);
    gr6(1) = Beta_1(1,6)/sigmasqalpha; gr6(2:nBeta) = Beta_1(2:nBeta,6)/tau_temp(6);
    gr7(1:nBeta) = Beta_2(1:nBeta,1)/tau_temp(7);
    gr8(1:nBeta) = Beta_2(1:nBeta,2)/tau_temp(8);
    gr9(1:nBeta) = Beta_2(1:nBeta,3)/tau_temp(9); 
    grs = [gr1;gr2;gr3;gr4;gr5;gr6;gr7;gr8;gr9];
    gr = gr + (nsubj-1)*grs;
end    