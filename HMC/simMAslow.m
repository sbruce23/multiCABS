function  TS0 = simMAslow(T0, nsubj, Sigma, cont, burn)


T=T0+burn;
d = size(Sigma,1);
TS0 = cell(nsubj,1);
for j=1:nsubj
    a1 = cont*(1-1.781*sin(pi*(j)/(2*nsubj)));
    a2 = cont*(1-1.781*cos(pi*(j)/(2*nsubj)));
    error = mvnrnd([0;0],Sigma,T+2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate starting terms
    a=[a1 -1; -1, a2];
    TS(:,1) = error(3,:)' + a*error(2,:)' + diag([0.5 -1.2])*error(1,:)';
    TS(:,2) = error(4,:)' + a*error(3,:)' + diag([0.5 -1.2])*error(2,:)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get the rest
    for t = 3:T
        TS(:,t) = error(t+2,:)' + a*error(t+1,:)' + diag([0.5 -1.2])*error(t,:)';
    end
    TS0{j} = TS';
 end    