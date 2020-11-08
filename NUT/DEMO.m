%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This file demonstrates the use of the program multiCABS(), which
%   implements the method discussed in 
%   "Conditional Adaptive Bayesian Spectral Analysis of Replicated Multivariate Time Series."
%
%   Content:
%   (1) Replicated bivariate slowly-varying time series  
%       (1a) get the bivariate slowly-varying time series
%       (1b) Set the options and hyperparameters for the sampler
%       (1c) Run the Bayesian estimation process
%       (1d) Obtain and plot the surface of spectra and coherence
%   (2) Replicated trivariate abrupt-changing time series 
%       (1a) get the trivariate abrupt-changing time series
%       (1b) Set the options and hyperparameters for the sampler
%       (1c) Run the Bayesian estimation process
%       (1d) Obtain and plot the surface of spectra and coherence

%% (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1a) Get the slowly-varying replicated time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load slow.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1b) Set the options and hyperparameters for the sampler
params = struct('nloop',2000, 'nwarmup',500, ...
                   'nexp_max',10, 'prob_mml',0.8 , 'nbasis',10, ...
                   'tau_up_limit',10^4, 'sigmasqalpha',10^4, 'init',3,...
                   'nfreq',100, 'verb',1, 'convdiag',0, 'ee', 0.1);   
                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1c) Run the Bayesian estimation process 
randn('state',19881205)
[spect_matrices, freq_hat, fit, ~] = multiCABS(zt,uu,params);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1d) Obtain and plot the estimated spectra and square coherence surfaces
[spect_slow, coh_slow] = multiCABS_surface(zt, spect_matrices, fit, params); 
bottom = 0; top = 7;
figure
subplot(1,3,1); 
meshc(linspace(0,1,N), freq_hat, real(spect_slow{1}));ax = gca; zlim([0,8])
ax.FontSize = 16; view(0,90);  title('$f_{11}$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
caxis manual
caxis([bottom top]);
colorbar('location','Manual','position',[0.35 0.12 0.01 0.8],'XTick',1:7); 
subplot(1,3,2); 
meshc(linspace(0,1,N), freq_hat, real(spect_slow{2})); 
ax = gca; zlim([0,8])
ax.FontSize = 16; view(0,90);  title('$f_{22}$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
caxis manual
caxis([bottom top]);
c=colorbar('location','Manual','position',[0.63 0.12 0.01 0.8],'XTick',1:7); 
subplot(1,3,3); 
meshc(linspace(0,1,N), freq_hat, real(coh_slow{1})); 
ax = gca; zlim([0,1])
ax.FontSize = 16; view(0,90);  title('$\rho_{21}^2$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
c=colorbar('location','Manual','position',[0.91 0.12 0.01 0.8]); 


%% (2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2a) Get the abrupt-change replicated time series 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load abrupt.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1b) Set the options and hyperparameters for the sampler
params = struct('nloop',2000, 'nwarmup',500, ...
                   'nexp_max',10, 'prob_mml',0.8 , 'nbasis',10, ...
                   'tau_up_limit',10^4, 'sigmasqalpha',10^4, 'init',3,...
                   'nfreq',50, 'verb',1, 'convdiag',0, 'ee', 0.1); 
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2c) Run the Bayesian estimation process 
randn('state',19881205)
[spect_matrices, freq_hat, fit, ~] = multiCABS(zt,uu,params);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2d) Obtain and plot the estimated spectra and square coherence surfaces
[spect, coh] = multiCABS_surface(zt, spect_matrices, fit, params); 

bottom = 0;
top  = 3;
figure
subplot(1,3,1); 
mesh(linspace(0,1,N), freq_hat,real(spect{1}));ax = gca; zlim([0,4])
ax.FontSize = 16; view(0,90);  title('$\hat{f}_{11}$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
%c=colorbar('location','Manual','position',[0.35 0.12 0.01 0.8],'XTick',[0, 1,1.5 2,2.5,3]); 
caxis manual
caxis([bottom top]);
subplot(1,3,2); 
mesh(linspace(0,1,N), freq_hat,real(spect{2})); 
ax = gca; zlim([0,4])
ax.FontSize = 16; view(0,90);  title('$\hat{f}_{22}$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
%c=colorbar('location','Manual','position',[0.63 0.12 0.01 0.8],'XTick',[0, 1,1.5 2,2.5,3]); 
caxis manual
caxis([bottom top]);
subplot(1,3,3); 
mesh(linspace(0,1,N), freq_hat,real(spect{3})); 
ax = gca; zlim([0,3])
ax.FontSize = 16; view(0,90);  title('$\hat{f}_{33}$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
caxis manual
caxis([bottom top]);
c=colorbar('location','Manual','position',[0.91 0.12 0.01 0.8]); 

bottom = 0;
top  = .5;
figure
subplot(1,3,1); 
mesh(linspace(0,1,N), freq_hat,(real(coh{1})));ax = gca; zlim([0,1])
ax.FontSize = 16; view(0,90);  title('$\hat{\rho}_{21}^2$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
%c=colorbar('location','Manual','position',[0.35 0.12 0.01 0.8],'XTick',[0, 0.2, 0.4, 0.6, 0.8, 1]); 
caxis manual
caxis([bottom top]);
subplot(1,3,2); 
mesh(linspace(0,1,N), freq_hat,(real(coh{2}))); 
ax = gca; zlim([0,1])
ax.FontSize = 16; view(0,90);  title('$\hat{\rho}_{31}^2$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
%c=colorbar('location','Manual','position',[0.63 0.12 0.01 0.8],'XTick',[0, 0.2, 0.4, 0.6, 0.8, 1]); 
caxis manual
caxis([bottom top]);
subplot(1,3,3); 
mesh(linspace(0,1,N), freq_hat,(real(coh{3}))); 
ax = gca; zlim([0,1])
ax.FontSize = 16; view(0,90);  title('$\hat{\rho}_{32}^2$','Interpreter','LaTex','Fontsize',24); xlabel('$\nu$','Interpreter','LaTex','Fontsize',24); ylabel('Freq','Fontsize',18); %zlim([0 8])
caxis manual
caxis([bottom top]);
c=colorbar('location','Manual','position',[0.91 0.12 0.01 0.8],'XTick',[0, 0.1, 0.2, 0.3, 0.4, 0.5]); 




