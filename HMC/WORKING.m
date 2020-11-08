
load svseed
randn('state',svseed(1))


ma11 = [0.6, 0 0; 0.2, -0.5, 0; 0.1 0.3 0.4];
ma12 = [0.3, 0 0; 0,  0.3 0 ; 0 0 0];
ma21 = [0.6, 0 0; 0.2, 0.5 0 ; -0.1 -0.3 0.4];
ma22 = [0.3, 0 0; 0, 0.3 0 ; 0 0 0];
sig  = [1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1];

T = 600;
N =  30;
cut = N/2;
zt = simMA(T,N,cut,ma11,ma12,ma21,ma22,sig,1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nloop = 2000;                   %number of total MCMC iterations
nwarmup = 1000;                 %number of warmup period
nexp_max = 5;                   %Input maximum number of segments

global dimen nobs nbasis nBeta sigmasqalpha tau_up_limit ...
            prob_mm1 umin M ee options    
        
if verLessThan('matlab','9.2')     
    options = optimset('Display','off','GradObj','on','Hessian','on','MaxIter',10000,...
    'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6);
else
    options = optimoptions(@fminunc,'Display','off','GradObj','on','Hessian','on','MaxIter',10000,...
        'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6,'Algorithm','trust-region');
end

nbasis = 10;                              %number of linear smoothing spline basis functions 
nBeta = nbasis + 1;                       %number of coefficients 
sigmasqalpha = 10^4;                      %smoothing parameter for alpha
tau_up_limit = 10^4;                      %the prior for smoothing parameters    
prob_mm1 = 0.8;                           %the probability of small jump, 1-prob of big jump
umin = 1;                                 %minimum number of observation in each segment
M = tau_up_limit;                         %scale value for momentum variable
ee = 0.1;                                %step size in Hamiltonian Monte Carlo
freq_hat=(1:floor((T-1)/2))'/(T);
nfreq_hat = 148;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% II) Run the estimation proceedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = size(zt{1});       
nobs = dim(1);
dimen = dim(2);
nsubj = size(zt,1);
ts = zt;
tausq = cell(nexp_max,1);       %big array for smoothing parameters
Beta = cell(nexp_max,1);        %big array for coefficients 
spect_hat = cell(nexp_max,1);
xi = cell(nexp_max,1);          %Cutpoint locations xi_1 is first cutpoint, xi_) is beginning of timeseries
nseg = cell(nexp_max,1);        %Number of observations in each segment
for j=1:nexp_max
    tausq{j}=ones(dimen^2,j,nloop+1);
    Beta{j} = zeros(nBeta,dimen^2,j,nloop+1);
    spect_hat{j} = zeros(dimen,dimen,nfreq_hat+1,j,nloop+1);
    xi{j}=ones(j,nloop+1);
	nseg{j}=ones(j,nloop+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilize the MCMC iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexp_curr = nexp_max*ones(nloop+1,1);     %big array for number of segment
nexp_curr(1) = 3;                         %initialize the number of segments

%initialize tausq
for j=1:nexp_curr(1)
	tausq{nexp_curr(1)}(:,j,1)=rand(dimen^2,1)*tau_up_limit;
end

%initilize the location of the changepoints for j=1:nexp_curr(1)
for j=1:nexp_curr(1)
    if nexp_curr(1)==1
        xi{nexp_curr(1)}(j,1) = nsubj;
        nseg{nexp_curr(1)}(j,1) = nsubj;
    else
        if j==1
			nposs = nsubj-nexp_curr(1)*umin+1;
			xi{nexp_curr(1)}(j,1) = umin + unidrnd(nposs)-1;
			nseg{nexp_curr(1)}(j,1) = xi{nexp_curr(1)}(j,1);
        elseif j>1 && j<nexp_curr(1)
			nposs=nsubj-xi{nexp_curr(1)}(j-1,1)-umin*(nexp_curr(1)-j+1)+1;
			xi{nexp_curr(1)}(j,1)=umin+unidrnd(nposs)+xi{nexp_curr(1)}(j-1,1)-1;
			nseg{nexp_curr(1)}(j,1)=xi{nexp_curr(1)}(j,1)-xi{nexp_curr(1)}(j-1,1);
        else
			xi{nexp_curr(1)}(j,1)=nsubj;
			nseg{nexp_curr(1)}(j,1)=xi{nexp_curr(1)}(j,1)-xi{nexp_curr(1)}(j-1,1);	
        end
    end
end

xi_temp = xi{nexp_curr(1)}(:,1);
tau_temp = tausq{nexp_curr(1)}(:,:,1);
for j=1:nexp_curr(1)
   	[Beta_mean, Beta_var,~] = postBeta(j, ts, tau_temp(:,j), xi_temp);
    if min(eig(0.5*(Beta_var+Beta_var')))<0
        Beta{nexp_curr(1)}(:,:,j,1) = reshape(mvnrnd(Beta_mean,0.5*(eye(nBeta*dimen^2))),nBeta,dimen^2);
    else    
        Beta{nexp_curr(1)}(:,:,j,1) = reshape(mvnrnd(Beta_mean,0.5*(Beta_var+Beta_var')),nBeta,dimen^2);
    end    
end

%jumping probabilities
epsilon=zeros(nloop,1);
met_rat=zeros(nloop,1);
bet_death = 0;
bet_birth = 0;
with = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:nloop
    tic;
    if(mod(p,50)==0)
       fprintf('iter: %g of %g \n' ,p, nloop)
       nexp_curr(p)
    end
    
    %========================
    %BETWEEN MODEL MOVE
    %========================
    kk = length(find(nseg{nexp_curr(p)}(:,p)>2*umin)); %Number of available segments
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Deciding on birth or death
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if kk==0 %Stay where you (if nexp_curr=1) or join segments if there are no available segments to cut
        if nexp_curr(p)==1
            nexp_prop = nexp_curr(p); %Stay 
            log_move_prop = 0;
            log_move_curr = 0;
        else
            nexp_prop = nexp_curr(p)-1; %death
            log_move_prop = 0;
            if nexp_prop==1
                log_move_curr = 1;
            else
                log_move_curr = log(0.5);
            end
        end
    else
        if nexp_curr(p)==1
            nexp_prop = nexp_curr(p) + 1; %birth
            log_move_prop = 0;
            if nexp_prop==nexp_max
                log_move_curr = 0;
            else
                log_move_curr = log(0.5);
            end
        elseif nexp_curr(p)==nexp_max
            nexp_prop = nexp_curr(p)-1;   %death
            log_move_prop = 0;
            if nexp_prop==1
                log_move_curr = 0;
            else
                log_move_curr = log(0.5);
            end
        else
            u = rand;
            if u<0.5
                nexp_prop = nexp_curr(p)+1; %birth
                if nexp_prop==nexp_max
                    log_move_curr = 0;
                    log_move_prop = log(0.5);
                else
                    log_move_curr=log(0.5);
                    log_move_prop=log(0.5);
                end
            else
                nexp_prop = nexp_curr(p)-1; %death
                if nexp_prop==1
                    log_move_curr = 0;
                    log_move_prop = log(0.5);
                else
                    log_move_curr = log(0.5);
                    log_move_prop = log(0.5);
                end
            end
        end
    end
    xi_curr_temp = xi{nexp_curr(p)}(:,p);
    Beta_curr_temp = Beta{nexp_curr(p)}(:,:,:,p);
    nseg_curr_temp = nseg{nexp_curr(p)}(:,p);
    tau_curr_temp = tausq{nexp_curr(p)}(:,:,p);

    if nexp_prop<nexp_curr(p)
        %Death step
        [met_rat(p),nseg_prop,xi_prop,tausq_prop,Beta_prop]= death(ts,nexp_curr(p),nexp_prop,...
		tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,log_move_curr,log_move_prop);
    elseif nexp_prop>nexp_curr(p)
        %Birth step
        [met_rat(p),nseg_prop,xi_prop,tausq_prop,Beta_prop]= birth(ts,nexp_curr(p),nexp_prop,...
		tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,log_move_curr,log_move_prop);
    else
        xi_prop=xi{nexp_curr(p)}(:,p);
        nseg_prop=nseg{nexp_curr(p)}(:,p);
        tausq_prop=tausq{nexp_curr(p)}(:,:,p);
        Beta_prop=Beta{nexp_curr(p)}(:,:,p);
        met_rat(p) = 1;
    end
    u = rand;
    if u<met_rat(p)
        if nexp_prop<nexp_curr(p)
            bet_death = bet_death + 1;
        elseif nexp_prop>nexp_curr(p)
            bet_birth = bet_birth + 1;
        end    
		nexp_curr(p+1)=nexp_prop;
		xi{nexp_curr(p+1)}(:,p+1)=xi_prop;
		nseg{nexp_curr(p+1)}(:,p+1)=nseg_prop;
		tausq{nexp_curr(p+1)}(:,:,p+1)=tausq_prop;
		Beta{nexp_curr(p+1)}(:,:,:,p+1)=Beta_prop;
    else
		nexp_curr(p+1)=nexp_curr(p);
		xi{nexp_curr(p+1)}(:,p+1)=xi{nexp_curr(p+1)}(:,p);
		nseg{nexp_curr(p+1)}(:,p+1)=nseg{nexp_curr(p+1)}(:,p);
        Beta{nexp_curr(p+1)}(:,:,:,p+1)=Beta{nexp_curr(p+1)}(:,:,:,p);
		tausq{nexp_curr(p+1)}(:,:,p+1)=tausq{nexp_curr(p+1)}(:,:,p);
    end
    
    %========================
    %WITHIN MODEL MOVE
    %========================
    
	%Drawing a new cut point and Beta simultaneously
    %update coeffiecient of linear basis function
    xi_curr_temp=xi{nexp_curr(p+1)}(:,p+1);
    Beta_curr_temp=Beta{nexp_curr(p+1)}(:,:,:,p+1);
    tau_temp=tausq{nexp_curr(p+1)}(:,:,p+1);
    nseg_curr_temp=nseg{nexp_curr(p+1)}(:,p+1);
    [epsilon(p),nseg_new,xi_prop,tausq_prop,Beta_prop,seg_temp]= ...
    within(ts,nexp_curr(p+1),tau_temp, xi_curr_temp,...
                                                  nseg_curr_temp, Beta_curr_temp);
    u = rand;
    if (u<epsilon(p)|| p==1) 
        with = with + 1;
        if nexp_curr(p+1)>1
            for j=seg_temp:seg_temp+1
                Beta{nexp_curr(p+1)}(:,:,j,p+1)=Beta_prop(:,:,j);
                xi{nexp_curr(p+1)}(j,p+1)=xi_prop(j);
                nseg{nexp_curr(p+1)}(j,p+1)=nseg_new(j);
            end
        else
            Beta{nexp_curr(p+1)}(:,:,p+1)=Beta_prop;
        end
    else
        Beta{nexp_curr(p+1)}(:,:,:,p+1)=Beta_curr_temp;
        xi{nexp_curr(p+1)}(:,p+1)=xi_curr_temp;
        nseg{nexp_curr(p+1)}(:,p+1)=nseg_curr_temp;
    end
    
    %Drawing tau
    for j=1:nexp_curr(p+1)
        for i=1:dimen^2
            if ismember(i,1:dimen + dimen*(dimen-1)/2)
                tau_a = nbasis/2;
                tau_b = sum(Beta{nexp_curr(p+1)}(2:end,i,j,p+1).^2)/2;
                u=rand;
                const1 = gamcdf(1/tau_up_limit,tau_a,1/tau_b);
                const2 = 1-u*(1-const1);
                tausq{nexp_curr(p+1)}(i,j,p+1) = 1/gaminv(const2,tau_a,1/tau_b);
            else   
                tau_a = nBeta/2;
                tau_b = sum(Beta{nexp_curr(p+1)}(1:nBeta,i,j,p+1).^2)/2;
                u=rand;
                const1 = gamcdf(1/tau_up_limit,tau_a,1/tau_b);
                const2 = 1-u*(1-const1);
                tausq{nexp_curr(p+1)}(i,j,p+1) = 1/gaminv(const2,tau_a,1/tau_b);
            end
        end     
    end   
    %==================================
    %Estimating Spectral Density
    %==================================
    [xx_r, xx_i] = lin_basis_func(freq_hat); %produce linear basis functions
    
    for j=1:nexp_curr(p+1)
        %getting the coefficients of linear basis functions
        g1 = Beta{nexp_curr(p+1)}(:,1:(dimen + dimen*(dimen-1)/2),j,p+1);
        g2 = Beta{nexp_curr(p+1)}(1:nBeta,(dimen + dimen*(dimen-1)/2 + 1):end,j,p+1);
        
        theta = zeros(dimen*(dimen-1)/2,(nfreq_hat+1));
        for i=1:dimen*(dimen-1)/2
            theta_real = xx_r * g1(:,i+dimen);
            theta_imag = xx_i * g2(:,i);
            theta(i,:) = (theta_real + sqrt(-1)*theta_imag).';
        end  
        
        delta_sq_hat = zeros(2,nfreq_hat+1);
        for q=1:dimen
            delta_sq_hat(q,:) = exp(xx_r * g1(:,q)).';
        end
        
        %produce the spectral density matrix
        for k=1:(nfreq_hat+1)
            TT = eye(dimen);
            TT(2,1) = -theta(1,k);
            if dimen==3
                TT(3,1) = -theta(2,k);
                TT(3,2) = -theta(3,k);
            end    
            spect_hat{nexp_curr(p+1)}(:,:,k,j,p+1) = ...
            inv(TT)*diag(delta_sq_hat(:,k))*inv(TT');  
        end
    end   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%true spectral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfreq_hat = 149;
%freq_hat=freq;
dim = size(zt{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f11 = zeros(nfreq_hat,N);
f22 = zeros(nfreq_hat,N);
f33 = zeros(nfreq_hat,N);
f21_real = zeros(nfreq_hat,N);
f21_imag = zeros(nfreq_hat,N);
f31_real = zeros(nfreq_hat,N);
f31_imag = zeros(nfreq_hat,N);
f32_real = zeros(nfreq_hat,N);
f32_imag = zeros(nfreq_hat,N);
f21 = zeros(nfreq_hat,N);
f31 = zeros(nfreq_hat,N);
f32 = zeros(nfreq_hat,N);
coh21_true=zeros(nfreq_hat,N);
coh31_true=zeros(nfreq_hat,N);
coh32_true=zeros(nfreq_hat,N);
tt=linspace(1,nobs,nobs)';
for i=1:N
    if i <= N/2
        beta =[ma11,ma12];
        spect_true=MAspectriv(beta,sig,freq_hat);
        f11(:,i) = squeeze(real(spect_true(1,1,:)));
        f22(:,i) = squeeze(real(spect_true(2,2,:)));
        f33(:,i) = squeeze(real(spect_true(3,3,:)));
        f21_real(:,i) = squeeze(real(spect_true(2,1,:)));
        f21_imag(:,i) = squeeze(imag(spect_true(2,1,:)));
        f31_real(:,i) = squeeze(real(spect_true(3,1,:)));
        f31_imag(:,i) = squeeze(imag(spect_true(3,1,:)));        
        f32_real(:,i) = squeeze(real(spect_true(3,2,:)));
        f32_imag(:,i) = squeeze(imag(spect_true(3,2,:)));              
        f21(:,i) = f21_real(:,i) + sqrt(-1)*f21_imag(:,i);
        f31(:,i) = f31_real(:,i) + sqrt(-1)*f31_imag(:,i);
        f32(:,i) = f32_real(:,i) + sqrt(-1)*f32_imag(:,i);
        coh21_true(:,i)=abs(f21(:,i)).^2./(f11(:,i).*f22(:,i));
        coh31_true(:,i)=abs(f31(:,i)).^2./(f11(:,i).*f33(:,i));
        coh32_true(:,i)=abs(f32(:,i)).^2./(f22(:,i).*f33(:,i));
    else
        beta =[ma21,ma22];
        spect_true=MAspectriv(beta,sig,freq_hat);
        f11(:,i) = squeeze(real(spect_true(1,1,:)));
        f22(:,i) = squeeze(real(spect_true(2,2,:)));
        f33(:,i) = squeeze(real(spect_true(3,3,:)));
        f21_real(:,i) = squeeze(real(spect_true(2,1,:)));
        f21_imag(:,i) = squeeze(imag(spect_true(2,1,:)));
        f31_real(:,i) = squeeze(real(spect_true(3,1,:)));
        f31_imag(:,i) = squeeze(imag(spect_true(3,1,:)));   
        f32_real(:,i) = squeeze(real(spect_true(3,2,:)));
        f32_imag(:,i) = squeeze(imag(spect_true(3,2,:)));           
        f21(:,i) = f21_real(:,i) + sqrt(-1)*f21_imag(:,i);
        f31(:,i) = f31_real(:,i) + sqrt(-1)*f31_imag(:,i);    
        f32(:,i) = f32_real(:,i) + sqrt(-1)*f32_imag(:,i);        
        coh21_true(:,i)=abs(f21(:,i)).^2./(f11(:,i).*f22(:,i)); 
        coh31_true(:,i)=abs(f31(:,i)).^2./(f11(:,i).*f33(:,i));     
        coh32_true(:,i)=abs(f32(:,i)).^2./(f22(:,i).*f33(:,i));
    end
end

s_11=zeros(nfreq_hat,N);
s_22=zeros(nfreq_hat,N);
s_33=zeros(nfreq_hat,N);
coh_21=zeros(nfreq_hat,N);
coh_31=zeros(nfreq_hat,N);
coh_32=zeros(nfreq_hat,N);

for p=1:nloop
    if(p>nwarmup)
        xi_curr=xi{nexp_curr(p)}(:,p);
        spec_hat_curr_11=squeeze(spect_hat{nexp_curr(p)}(1,1,:,:,p));
        spec_hat_curr_22=squeeze(spect_hat{nexp_curr(p)}(2,2,:,:,p));
        spec_hat_curr_33=squeeze(spect_hat{nexp_curr(p)}(3,3,:,:,p));
        spec_hat_curr_21=abs(squeeze(spect_hat{nexp_curr(p)}(2,1,:,:,p))).^2./...
            (squeeze(spect_hat{nexp_curr(p)}(1,1,:,:,p)).*squeeze(spect_hat{nexp_curr(p)}(2,2,:,:,p)));
        spec_hat_curr_31=abs(squeeze(spect_hat{nexp_curr(p)}(3,1,:,:,p))).^2./...
            (squeeze(spect_hat{nexp_curr(p)}(1,1,:,:,p)).*squeeze(spect_hat{nexp_curr(p)}(3,3,:,:,p)));
        spec_hat_curr_32=abs(squeeze(spect_hat{nexp_curr(p)}(3,2,:,:,p))).^2./...
            (squeeze(spect_hat{nexp_curr(p)}(2,2,:,:,p)).*squeeze(spect_hat{nexp_curr(p)}(3,3,:,:,p)));
        for j=1:nexp_curr(p)
            if(j==1)
                s_11(:,1:xi_curr)=s_11(:,1:xi_curr(j))+repmat(spec_hat_curr_11(:,j),1,xi_curr(j))/(nloop-nwarmup);
                s_22(:,1:xi_curr)=s_22(:,1:xi_curr(j))+repmat(spec_hat_curr_22(:,j),1,xi_curr(j))/(nloop-nwarmup);
                s_33(:,1:xi_curr)=s_33(:,1:xi_curr(j))+repmat(spec_hat_curr_33(:,j),1,xi_curr(j))/(nloop-nwarmup);
                coh_21(:,1:xi_curr)=coh_21(:,1:xi_curr(j))+repmat(spec_hat_curr_21(:,j),1,xi_curr(j))/(nloop-nwarmup);
                coh_31(:,1:xi_curr)=coh_31(:,1:xi_curr(j))+repmat(spec_hat_curr_31(:,j),1,xi_curr(j))/(nloop-nwarmup);
                coh_32(:,1:xi_curr)=coh_32(:,1:xi_curr(j))+repmat(spec_hat_curr_32(:,j),1,xi_curr(j))/(nloop-nwarmup);

            else
                s_11(:,xi_curr(j-1)+1:xi_curr(j))=s_11(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_11(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
                s_22(:,xi_curr(j-1)+1:xi_curr(j))=s_22(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_22(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
                s_33(:,xi_curr(j-1)+1:xi_curr(j))=s_33(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_33(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
                coh_21(:,xi_curr(j-1)+1:xi_curr(j))=coh_21(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_21(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
                coh_31(:,xi_curr(j-1)+1:xi_curr(j))=coh_31(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_31(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
                coh_32(:,xi_curr(j-1)+1:xi_curr(j))=coh_32(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_32(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);                    
            end
        end
    end
end