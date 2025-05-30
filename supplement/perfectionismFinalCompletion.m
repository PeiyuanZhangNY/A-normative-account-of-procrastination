%% when beta>lambda, final completion 0 or 1; 
NMonte = 10000;
deltas = 0.01;k=1;Tmax = 10;Tmin=2;TAll = randi([Tmin Tmax],1,NMonte);
gammaAll = rand(1,NMonte);
lambdaAll = 10*rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
betaAll = lambdaAll+(10-lambdaAll).*rand(1,NMonte); % beta>lambda
%betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
%betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda

alphaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll(1,NMonte).*rand(1,NMonte);

finalprop = nan(1,NMonte);
for nMonte =1:NMonte
    nMonte
    [temp1,] = OptActStateSeq(deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 

    finalprop(1,nMonte)=sum(temp1);
end

histogram(finalprop)
%save('perfectionism_betalargerlambda.mat')

%%%%%%%%%% conclusion: finalprop is either 0 or 1. 

%% when beta<lambda, final completion is between 0 and 1
NMonte = 10000;
deltas = 0.01;k=1;Tmax = 10;Tmin=2;TAll = randi([Tmin Tmax],1,NMonte);
gammaAll = rand(1,NMonte);
lambdaAll = 10*rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
%betaAll = lambdaAll+(10-lambdaAll).*rand(1,NMonte); % beta>lambda
betaAlllog = log10(0.1)+(log10(lambdaAll)-log10(0.1)).*rand(1,NMonte);
betaAll = 10.^betaAlllog; % beta <lambda

alphaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll(1,NMonte).*rand(1,NMonte);

finalprop = nan(1,NMonte);
for nMonte =1:NMonte
    nMonte
    [temp1,] = OptActStateSeq(deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 

    finalprop(1,nMonte)=sum(temp1);
end

histogram(finalprop)
%save('perfectionism_betasmallerlambda.mat')
%%%%%%%%%% conclusion: finalprop is between 0 and 1,[0,1]. 