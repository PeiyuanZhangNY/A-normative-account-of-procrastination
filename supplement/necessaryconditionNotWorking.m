% necessary condition to find not working at all. 

NMonte = 10000;
deltas = 0.01;k=1;Tmax = 10;Tmin=2;
TAll = randi([Tmin Tmax],1,NMonte);
gammaAll = rand(1,NMonte);
lambdaAll = 10*rand(1,NMonte); 
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda

%% when J<alpha-c_1, J>0, so alpha-c_1>0 
alphaAll = rand(1,NMonte);
c1All = alphaAll.*rand(1,NMonte);
JAll = (alphaAll-c1All).*rand(1,NMonte);

finalprop1 = nan(1,NMonte);

for nMonte =1:NMonte
    nMonte
    [temp,] = OptActStateSeq(deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
    finalprop1(nMonte)=sum(temp);
end

%%%%% the result: when J<alpha-c_1, all is to work on the task (finalprop>0), in some situations 
%%%%% partially done, some situations completion;

%% when J>alpha-c_1, J>0
alphaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
JAll = max((alphaAll-c1All)+(5-(alphaAll-c1All)).*rand(1,NMonte),0);

finalprop2 = nan(1,NMonte);

for nMonte =1:NMonte
    nMonte
    [temp,] = OptActStateSeq(deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
    finalprop2(nMonte)=sum(temp);
end
%%%%%%% when J>alpha-c_1, final proportion completed <=1 (existing of 0's). So
%%%%%%% J>alpha-c_1 is a necessary condition to get final proportion
%%%%%%% completed as 0.



