function [spect_hat, freq_hat, diagparams, fitparams] = multiCABS(X,uu,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function for the Conditional Adaptive Bayesian Spectral Analysis of
% Replicated Multivariate Time Series
% Does the MCMC iterations for bivariate or trivariate time series

%   Input:
%       1) zt - is a cell contains N replicated P dimensional time series
%       2) uu - a N-dimensional vector contains the covariate values 
%       3) varargin - tuning parameters, which can be set by the program OptsMultiCABS.m  
%       See the documentation in that program for details and default values. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  I) Extract information from the option parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If params is empty, then use the default paraemters
if  nargin==1
    params = OptsMultiSpect();
else
    params = varargin{1};
end

nloop = params.nloop;                       %number of total MCMC iterations
nwarmup = params.nwarmup;                   %number of warmup period
nexp_max = params.nexp_max;                 %Input maximum number of segments

global dimen nobs nbasis nBeta sigmasqalpha tau_up_limit ...
            prob_mm1 umin M options  
        
if verLessThan('matlab','9.2')     
    options = optimset('Display','off','GradObj','on','Hessian','on','MaxIter',10000,...
    'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6);
else
    options = optimoptions(@fminunc,'Display','off','GradObj','on','Hessian','on','MaxIter',10000,...
        'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6,'Algorithm','trust-region');
end

nbasis = params.nbasis;                   %number of linear smoothing spline basis functions 
nBeta = nbasis + 1;                       %number of coefficients 
sigmasqalpha = params.sigmasqalpha;       %smoothing parameter for alpha
tau_up_limit = params.tau_up_limit;       %the prior for smoothing parameters    
prob_mm1 = params.prob_mml;               %the probability of small jump, 1-prob of big jump
umin = 1;                                 %minimum number of observation in each segment
M = tau_up_limit;                         %scale value for momentum variable
nfreq_hat = params.nfreq;
freq_hat=(0:nfreq_hat)'/(2*nfreq_hat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% II) Run the estimation proceedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = size(X{1});       
nobs = dim(1);
dimen = dim(2);
nsubj = size(X,1);
[uu,I] = sort(uu); 
ts = X(I);
u_unique = unique(uu);
for i=1:length(u_unique)
    uu_cut(i) = max(find(uu==u_unique(i)));
end
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
nexp_curr(1) = params.init;               %initialize the number of segments

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
   	[Beta_mean, Beta_var,~] = postBeta(j, ts, tau_temp(:,j), xi_temp, uu_cut);
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
		tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,log_move_curr,log_move_prop, uu_cut);
    elseif nexp_prop>nexp_curr(p)
        %Birth step
        [met_rat(p),nseg_prop,xi_prop,tausq_prop,Beta_prop]= birth(ts,nexp_curr(p),nexp_prop,...
		tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,log_move_curr,log_move_prop, uu_cut);
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
                                               nseg_curr_temp, Beta_curr_temp,uu_cut);                                      
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
    
    tms(p) = toc;
    if params.verb ==1
        fprintf('sec / min hr: %g %g %g \n' ,[tms(p),sum(tms(1:p))/60, sum(tms(1:p))/(60*60)]')          
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  III) Collect outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fitparams = struct('nloop', nloop, ...
                   'nwarmup', nwarmup, ...
                   'timeMean', mean(tms), ...
                   'timeMax', max(tms(2:end)), ...
                   'timeMin', min(tms),...
                   'timeStd', std(tms(2:end)));
diagparams = struct('tausq', tausq,...
                    'Beta', Beta,...
                    'xi', xi,...
                    'nseg', nseg,...
                    'nexp_curr', nexp_curr,... 
                    'epsilon', epsilon,...                 
                    'bet_birth', bet_birth,...
                    'bet_death', bet_death,...
                    'with', with,...
                    'time', tms); 
end
