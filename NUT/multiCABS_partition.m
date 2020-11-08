function[post_partitions] = multiCABS_partition(zt, fit, params) 

dim = size(zt);nobs = dim(1);
post_partitions = (histc(fit(1).nexp_curr(params.nwarmup+1:params.nloop),1:params.nexp_max)'/length(params.nwarmup+1:params.nloop))';

figure
histogram(fit(1).nexp_curr(params.nwarmup+1:params.nloop),'FaceAlpha',1,...
    'Normalization','probability','BinMethod','integers','BinWidth',.5,'BinLimits',[1,10])
title('Histogram of the Number of Partitions')
xlabel('The number of partitions')
ylabel('Probability')
for j=1:params.nexp_max
	kk=find(fit(1).nexp_curr(params.nwarmup+1:params.nloop)==j);
    if ~isempty(kk) && j>1
        figure
        hold
        title(['Plot of Partition Points Given ',int2str(j), ' Segments'])
        for k=1:j-1
            plot(fit(j).xi(k,kk+params.nwarmup))
        end
        for k=1:j-1
            figure
            hold
            title(['Histogram of Location of Partition ', int2str(k), ',' ' Given ',int2str(j), ' Segments'])
            histogram(fit(j).xi(k,kk+params.nwarmup),'FaceAlpha',1,...
                'Normalization','probability','BinMethod','integers','BinWidth',20,'BinLimits',[1,nobs])
        end
    end
end
