function  TS0 = simMA(T0, nsubj, cut, ma11, ma12, ma21, ma22, Sigma, burn)


T=T0+burn;
d = size(Sigma,1);
TS0 = cell(nsubj,1);
for j=1:nsubj
    if j<=cut
        error = mvnrnd(zeros(d,1),Sigma,T+2);
        TSS(:,1) = error(3,:) + (ma11*error(2,:)')' + (ma12*error(1,:)')';
        TSS(:,2) = error(4,:) + (ma11*error(3,:)')' + (ma12*error(2,:)')';
        for t = 3:T
            TSS(:,t) = error(t+2,:) + (ma11*error(t+1,:)')' + (ma12*error(t,:)')';
        end
        TSS0 = TSS(:,(burn+1):end);
        TS0{j} = TSS0';
    else
        error = mvnrnd(zeros(d,1),Sigma,T+2);
        TSS(:,1) = error(3,:) + (ma21*error(2,:)')' + (ma22*error(1,:)')';
        TSS(:,2) = error(4,:) + (ma21*error(3,:)')' + (ma22*error(2,:)')';
        for t = 3:T
            TSS(:,t) = error(t+2,:) + (ma21*error(t+1,:)')' + (ma22*error(t,:)')';
        end 
        TSS0 = TSS(:,(burn+1):end);
        TS0{j} = TSS0';
    end    
 end    

