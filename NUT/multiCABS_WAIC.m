function[waic_1, waic_2] = multiCABS_WAIC(X, fit, params) 


dim = size(X{1});       
nobs = dim(1);
nsubj = size(X,1);

lppd_tmp = zeros(floor(nobs/2)+1, nsubj, params.nloop-params.nwarmup);

for p=1:params.nloop
    if p>params.nwarmup
        xi_curr=fit(fit(1).nexp_curr(p)).xi(:,p);
        for j=1:fit(1).nexp_curr(p)
            if(j==1)
                for k=1:xi_curr(j)
                    lppd_tmp(:,k,p-params.nwarmup) = lppd(X{k},fit(fit(1).nexp_curr(p)).Beta(:,:,j,p)); 
                end
            else
                for k=xi_curr(j-1)+1:xi_curr(j)
                    lppd_tmp(:,k,p-params.nwarmup) = lppd(X{k},fit(fit(1).nexp_curr(p)).Beta(:,:,j,p))'; 
                end
            end    
        end   
    end
end
lppd_value = sum(sum(log(mean(lppd_tmp,3))));
pwaic_1 = sum(sum(log(mean(lppd_tmp,3)) - mean(log(lppd_tmp),3)));
pwaic_2 = sum(sum(var(log(lppd_tmp),0,3)));
waic_1 =  lppd_value -  pwaic_1;
waic_2 =  lppd_value -  pwaic_2;