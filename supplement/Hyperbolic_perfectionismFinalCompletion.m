%% when beta>lambda, final completion 0 or 1; 
NMonte = 10000;
k=1;Tmax = 10;Tmin=2;TAll = randi([Tmin Tmax],1,NMonte);
lambdaAll = 10*rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
betaAll = lambdaAll+(10-lambdaAll).*rand(1,NMonte); % beta>lambda
%betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
%betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda

alphaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll(1,NMonte).*rand(1,NMonte);

logk = -9+(1-(-9))*rand(1,NMonte);
k_DR = exp(logk);

finalprop = nan(1,NMonte);
for nMonte =1:NMonte
    nMonte
    [temp1,] = HyperbolicDiscountingFun(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(1,nMonte),betaAll(nMonte),k_DR(nMonte)); 

    finalprop(1,nMonte)=sum(temp1);
end

hist(finalprop)
%save('hyperbolic_perfectionism_betalargerlambda.mat')
sum(finalprop > 0 & finalprop < 1) % 1 in 1000 simulations (it returns 9 cases, but 8 of them are 1's)

%%%%%%%%% conclusion: finalprop is mostly 0 or 1, with a few exceptions. 

%% when beta<lambda, final completion is between 0 and 1
NMonte = 10000;
k=1;Tmax = 10;Tmin=2;TAll = randi([Tmin Tmax],1,NMonte);
lambdaAll = 10*rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
%betaAll = lambdaAll+(10-lambdaAll).*rand(1,NMonte); % beta>lambda
betaAlllog = log10(0.1)+(log10(lambdaAll)-log10(0.1)).*rand(1,NMonte);
betaAll = 10.^betaAlllog; % beta <lambda

alphaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll(1,NMonte).*rand(1,NMonte);

logk = -9+(1-(-9))*rand(1,NMonte);
k_DR = exp(logk);

finalprop = nan(1,NMonte);
for nMonte =1:NMonte
    nMonte
    [temp1,] = HyperbolicDiscountingFun(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(1,nMonte),betaAll(nMonte),k_DR(nMonte)); 

    finalprop(1,nMonte)=sum(temp1);
end

histogram(finalprop)
%save('hyperbolic_perfectionism_betasmallerlambda.mat')
sum(finalprop > 0 & finalprop < 1) % 263-6 in 1000 simulations (6 cases are 1's). % sort(finalprop(finalprop > 0 & finalprop < 1));
%%%%%%%%%% conclusion: finalprop is between 0 and 1,[0,1]. 