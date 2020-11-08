function [f,gr,h] = Beta_derive_replicate(x, xobs, tau_temp, sigmasqalpha, nbasis)

global dimen
nsubj = size(xobs,1);
f = 0;
gr = 0;
h = 0;
for j=1:nsubj
   [f_tmp, gr_tmp, h_tmp] = Beta_derive_single(x, xobs{j}, tau_temp, sigmasqalpha, nbasis); 
   f = f + f_tmp;
   gr = gr + gr_tmp;
   h = h + h_tmp;
end  

%initilize Beta_1 and Beta_2: Beta_2 is for imaginary components
nBeta = nbasis + 1;
Beta_1 = zeros(nBeta,(dimen + dimen*(dimen-1)/2));
Beta_2 = zeros(nBeta,dimen*(dimen-1)/2);

x = reshape(x,nBeta,dimen^2);
Beta_temp = x;
Beta_1(:,:) = Beta_temp(:,1:(dimen + dimen*(dimen-1)/2));
Beta_2(:,:) = Beta_temp(1:nBeta,(dimen + dimen*(dimen-1)/2 + 1): end);

if dimen==2


    f = f + (nsubj-1)*(0.5.*(Beta_1(1,1)* Beta_1(1,1)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,1)'*Beta_1(2:nBeta,1))/tau_temp(1)) -...
            (0.5.*(Beta_1(1,2)* Beta_1(1,2)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,2)'*Beta_1(2:nBeta,2))/tau_temp(2)) -...
            (0.5.*(Beta_1(1,3)* Beta_1(1,3)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,3)'*Beta_1(2:nBeta,3))/tau_temp(3)) -...
            (0.5.*(Beta_2(1:nBeta,1)'*Beta_2(1:nBeta,1))/tau_temp(4));
       
    gr1 = zeros(nBeta,1); gr2 = zeros(nBeta,1); gr3 = zeros(nBeta,1); gr4 = zeros(nBeta,1);
    gr1(1) = Beta_1(1,1)/sigmasqalpha; gr1(2:nBeta,1) = Beta_1(2:nBeta,1)/tau_temp(1);
    gr2(1) = Beta_1(1,2)/sigmasqalpha; gr2(2:nBeta,1) = Beta_1(2:nBeta,2)/tau_temp(2);
    gr3(1) = Beta_1(1,3)/sigmasqalpha; gr3(2:nBeta,1) = Beta_1(2:nBeta,3)/tau_temp(3);
    gr4(1:nBeta,1) = Beta_2(1:nBeta,1)/tau_temp(4);   
    grs = [gr1;gr2;gr3;gr4];
    gr = gr - (nsubj-1)*grs;
    
    
    h11(1,1)=1/sigmasqalpha; h11(2:nBeta,2:nBeta)= 1/tau_temp(1)*eye(nbasis);
    h22(1,1)=1/sigmasqalpha; h22(2:nBeta,2:nBeta)= 1/tau_temp(2)*eye(nbasis);
    h33(1,1)=1/sigmasqalpha; h33(2:nBeta,2:nBeta)= 1/tau_temp(3)*eye(nbasis);
    h44(1:nBeta,1:nBeta)=1/tau_temp(4)*eye(nBeta);
    h42 = zeros(nBeta,nBeta); h32 = zeros(nBeta,nBeta);
    h24=h42'; h23=h32';
    h1 = [h11,zeros(nBeta,2*nBeta+nBeta)];
    h2 = [zeros(nBeta,nBeta),h22,-h23,-h24];
    h3 = [zeros(nBeta,nBeta),-h32,h33,zeros(nBeta,nBeta)];
    h4 = [zeros(nBeta,nBeta),-h42,zeros(nBeta,nBeta),h44]; 
    hs = [h1;h2;h3;h4];
    h = h - (nsubj-1)*hs;
else
   f = f + (nsubj-1)*(0.5.*(Beta_1(1,1)* Beta_1(1,1)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,1)'*Beta_1(2:nBeta,1))/tau_temp(1)) -...
          (0.5.*(Beta_1(1,2)* Beta_1(1,2)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,2)'*Beta_1(2:nBeta,2))/tau_temp(2)) -...
          (0.5.*(Beta_1(1,3)* Beta_1(1,3)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,3)'*Beta_1(2:nBeta,3))/tau_temp(3)) -...
          (0.5.*(Beta_1(1,4)* Beta_1(1,4)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,4)'*Beta_1(2:nBeta,4))/tau_temp(4)) -...
          (0.5.*(Beta_1(1,5)* Beta_1(1,5)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,5)'*Beta_1(2:nBeta,5))/tau_temp(5)) -...
          (0.5.*(Beta_1(1,6)* Beta_1(1,6)')/sigmasqalpha + 0.5.*(Beta_1(2:nBeta,6)'*Beta_1(2:nBeta,6))/tau_temp(6)) -...
          (0.5.*(Beta_2(1:nBeta,1)'*Beta_2(1:nBeta,1))/tau_temp(7))-...
          (0.5.*(Beta_2(1:nBeta,2)'*Beta_2(1:nBeta,2))/tau_temp(8))-...
          (0.5.*(Beta_2(1:nBeta,3)'*Beta_2(1:nBeta,3))/tau_temp(9));
    
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
    gr = gr - (nsubj-1)*grs;
    
    h11(1,1)=1/sigmasqalpha; h11(2:nBeta,2:nBeta)=1/tau_temp(1)*eye(nbasis);
    h22(1,1)=1/sigmasqalpha; h22(2:nBeta,2:nBeta)= 1/tau_temp(2)*eye(nbasis);
    h33(1,1)=1/sigmasqalpha; h33(2:nBeta,2:nBeta)= 1/tau_temp(3)*eye(nbasis);
    h44(1,1)=1/sigmasqalpha; h44(2:nBeta,2:nBeta)=1/tau_temp(4)*eye(nbasis);
    h55(1,1)=1/sigmasqalpha; h55(2:nBeta,2:nBeta)= 1/tau_temp(5)*eye(nbasis);
    h66(1,1)=1/sigmasqalpha; h66(2:nBeta,2:nBeta)= 1/tau_temp(6)*eye(nbasis);
    h77(1:nBeta,1:nBeta)=1/tau_temp(7)*eye(nBeta);
    h88(1:nBeta,1:nBeta)=1/tau_temp(8)*eye(nBeta);
    h99(1:nBeta,1:nBeta)=1/tau_temp(9)*eye(nBeta);
    h42 = zeros(nBeta,nBeta);
    h53 = zeros(nBeta,nBeta);
    h56 = zeros(nBeta,nBeta);
    h59 = zeros(nBeta,nBeta);
    h63 = zeros(nBeta,nBeta); 
    h68 = zeros(nBeta,nBeta);
    h72 = zeros(nBeta,nBeta);
    h83 = zeros(nBeta,nBeta);
    h93 = zeros(nBeta,nBeta);
    h98 = zeros(nBeta,nBeta);
    h24=h42'; h35=h53'; h65=h56'; h95=h59'; h36=h63';
    h86=h68'; h27=h72'; h38=h83'; h39=h93'; h89=h98';
    ze = zeros(nBeta,nBeta);
    h1 = [h11,repmat(ze,1,8)];
    h2 = [ze,h22,ze,-h24,ze,ze,-h27,ze,ze];
    h3 = [ze,ze,h33,ze,-h35,-h36,ze,-h38,-h39];
    h4 = [ze,-h42,ze,h44,repmat(ze,1,5)];  
    h5 = [ze,ze,-h53,ze,h55,h56,ze,ze,h59];
    h6 = [ze,ze,-h63,ze,h65,h66,ze,h68,ze];
    h7 = [ze,-h72,ze,ze,ze,ze,h77,ze,ze];
    h8 = [ze,ze,-h83,ze,ze,h86,ze,h88,h89];
    h9 = [ze,ze,-h93,ze,h95,ze,ze,h98,h99];
    hs = [h1;h2;h3;h4;h5;h6;h7;h8;h9];
    h = h - (nsubj-1)*hs;
end    