function[A,nseg_new,xi_prop,tau_prop,Beta_prop,seg_temp]=...
        within(ts,nexp_temp,tau_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp, uu_cut)
global nobs dimen nBeta prob_mm1 umin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does the within-model move in the paper
%
%   Input:
%       1) ts - a cell contains N replicated P dimensional time series
%       2) nexp_temp - number of segments
%       3) tau_temp -  vector of smoothing parameters
%       4) xi_curr_temp - current partitions
%       5) nseg_curr_temp - current number of observations in each segment
%       6) Beta_curr_temp - current vector of coefficients
%       7) uu_cut - group partitions
%   Main Outputs:
%       1) PI - acceptance probability
%       2) nseg_new - new number of observations in each segment
%       3) xi_prop - proposed partitions
%       4) tau_prop - proposed smoothing parameters
%       5) Beta_prop - proposed coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xi_prop = xi_curr_temp;
Beta_prop = Beta_curr_temp;
nseg_new = nseg_curr_temp;
tau_prop = tau_temp;


if nexp_temp>1
    %*********************************************************
    % If contains more than one segments
    %*********************************************************

    seg_temp = unidrnd(nexp_temp-1);  %Drawing Segment to cut
    u = rand;
    cut_poss_curr = xi_curr_temp(seg_temp);
    nposs_prior = nseg_curr_temp(seg_temp) + nseg_curr_temp(seg_temp+1) - 2*umin+1;
    
    %Determing if the relocation is a big jump or small jump
    if u<prob_mm1
        if nseg_curr_temp(seg_temp)==umin && nseg_curr_temp(seg_temp+1)==umin
            nposs=1; %Number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp)- 1 + new_index;
        elseif nseg_curr_temp(seg_temp)==umin
			nposs=2; %Number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp)- 1 + new_index;
        elseif nseg_curr_temp(seg_temp+1)==umin
			nposs=2; %Number of possible locations for new cutpoint
			new_index = unidrnd(nposs); %Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp) + 1 - new_index;
        else
			nposs=3;% Number of possible locations for new cutpoint
			new_index = unidrnd(nposs);%Drawing index of new cutpoint 
			cut_poss_new = xi_curr_temp(seg_temp) - 2 + new_index;
        end
    else
		new_index=unidrnd(nposs_prior);%
		cut_poss_new = sum(nseg_curr_temp(1:seg_temp-1))- 1 + umin+new_index;
    end
    
    xi_prop(seg_temp)=cut_poss_new;
    if seg_temp>1
        %Number of observations in lower part of new cutpoin
		nseg_new(seg_temp) = xi_prop(seg_temp) - xi_curr_temp(seg_temp-1); 
    else
		nseg_new(seg_temp) = xi_prop(seg_temp);
    end
    %Number of observations in upper part of new cutpoint
	nseg_new(seg_temp+1) = nseg_curr_temp(seg_temp) + nseg_curr_temp(seg_temp+1) - nseg_new(seg_temp);
    
    %=========================================================================================
    %Evaluating the cut Proposal density for the cut-point at the cureent and proposed values
    %=========================================================================================
    if(abs(cut_poss_new-cut_poss_curr)>1)
        log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior);
        log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior);
    elseif nseg_curr_temp(seg_temp)==umin && nseg_curr_temp(seg_temp+1)==umin
        log_prop_cut_prop=0;
        log_prop_cut_curr=0;
    else
        if (nseg_curr_temp(seg_temp)==umin || nseg_curr_temp(seg_temp+1)==umin)
           %log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior)+log(1/2)+log(prob_mm1);
           log_prop_cut_prop=log(1/2)+log(prob_mm1);
        else
           %log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior)+log(1/3)+log(prob_mm1); 
           log_prop_cut_prop=log(1/3)+log(prob_mm1); 
        end
        if(nseg_new(seg_temp)==umin || nseg_new(seg_temp+1)==umin)
           %log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior)+log(1/2)+log(prob_mm1);
           log_prop_cut_curr=log(1/2)+log(prob_mm1);
        else
           %log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior)+log(1/3)+log(prob_mm1); 
           log_prop_cut_curr=log(1/3)+log(prob_mm1); 
        end
    end
    
    select = 1:dimen^2;

    %==========================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the current values
    %==========================================================================
    loglike_curr = 0;
    for j=seg_temp:seg_temp+1 
            if j>1
                yobs_tmp = ts((xi_curr_temp(j-1)+1):xi_curr_temp(j),:);
            else
                yobs_tmp = ts(1:xi_curr_temp(j),:);
            end
         %Loglikelihood at current values
         [log_curr_spec_dens] = whittle_like_replicate(yobs_tmp,Beta_curr_temp(:,:,j));  
         loglike_curr = loglike_curr + log_curr_spec_dens/(nobs/2);
    end
    
    %=================================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the proposed values
    %=================================================================================
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For coefficient of basis functions that are different across two segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loglike_prop = 0;
    yobs_tmp_2 = cell(2,1);
    for j=seg_temp:seg_temp+1

        [Beta_out, yobs_tmp] = NUT_sampler(j, ts, tau_prop(:,j), Beta_prop(:,:,j), xi_prop, uu_cut);                            
        if j==seg_temp
            yobs_tmp_2{1} = yobs_tmp;
        else
            yobs_tmp_2{2} = yobs_tmp;
        end                    
        
        Beta_prop(:,select,j) = reshape(Beta_out,nBeta,length(select));       
    end   
    
    %Loglikelihood at proposed values
    for j=seg_temp:seg_temp+1
        if j==seg_temp
            [log_curr_spec_dens] = whittle_like_replicate(yobs_tmp_2{1},Beta_prop(:,:,j)); 
        else
            [log_curr_spec_dens] = whittle_like_replicate(yobs_tmp_2{2},Beta_prop(:,:,j));
        end
        loglike_prop = loglike_prop+log_curr_spec_dens/(nobs/2);    
    end

    %proposal density
    log_proposal_curr =  log_prop_cut_curr;
    log_proposal_prop =  log_prop_cut_prop;
    
    %target density
    log_prior_cut_prop=0;
	log_prior_cut_curr=0;
    for k=1:nexp_temp-1
        if k==1
            log_prior_cut_prop=-log(nobs-(nexp_temp-k+1)*umin+1);
			log_prior_cut_curr=-log(nobs-(nexp_temp-k+1)*umin+1);
		else
			log_prior_cut_prop=log_prior_cut_prop - log(nobs-xi_prop(k-1)-(nexp_temp-k+1)*umin+1);
			log_prior_cut_curr=log_prior_cut_curr - log(nobs-xi_curr_temp(k-1)-(nexp_temp-k+1)*umin+1);
        end
    end
    log_target_prop = loglike_prop + log_prior_cut_prop;
    log_target_curr = loglike_curr + log_prior_cut_curr;
    
else    
    
    %*********************************************************
    % If contains only one segment
    %*********************************************************
    nseg_new = nobs;
    seg_temp = 1;
    %=================================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the proposed values
    %=================================================================================
    
    [Beta_out,~] = NUT_sampler(1, ts, tau_temp, Beta_curr_temp, xi_curr_temp, uu_cut);
    Beta_prop(:,:,1) = reshape(Beta_out,nBeta,dimen^2);
  

    %Loglike at proposed values
    [loglike_prop] = whittle_like_replicate(ts,Beta_prop(:,:,1));
    
    %Loglike at current values
    [loglike_curr] = whittle_like_replicate(ts,Beta_curr_temp(:,:,1));
    
    log_target_prop = loglike_prop/(nobs/2);
    log_target_curr = loglike_curr/(nobs/2);
    log_proposal_curr =  0;
    log_proposal_prop =  0;
end
%loglike_prop-loglike_curr
%*************************************************************
%Calculations acceptance probability	
%*************************************************************
A = min(1,exp(log_target_prop-log_target_curr+log_proposal_curr-log_proposal_prop));
