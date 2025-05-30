
deltas = 0.01; k=1;u=0;
eta=0.8; gamma = 0.2;
alpha=1; beta=1;
deltalambda=0.4;
lambdaVec = 0.2:deltalambda:4;
T = 7;

WithoptActSeqMatrix=nan(T,length(lambdaVec));
U_task_with=nan(1,length(lambdaVec));
cost_with=nan(1,length(lambdaVec));
Utotal_with=nan(1,length(lambdaVec));
meanWithoptActSeqMatrix = nan(1,length(lambdaVec));

for lambdaIdx=1:length(lambdaVec)
    [WithoptActSeqMatrix(:,lambdaIdx), ]=OptActStateSeq( deltas,T,k,u,eta,lambdaVec(lambdaIdx),gamma,alpha,beta);    
end

for lambdaIdx=1:length(lambdaVec)  
    meanWithoptActSeqMatrix(lambdaIdx) = nanmean(WithoptActSeqMatrix(:,lambdaIdx));
    U_task_with(lambdaIdx)=alpha*(k*nansum(WithoptActSeqMatrix(:,lambdaIdx)))^beta;
    cost_with(lambdaIdx)=nansum(eta*WithoptActSeqMatrix(:,lambdaIdx).^lambdaVec(lambdaIdx));
    Utotal_with(lambdaIdx) = U_task_with(lambdaIdx)-cost_with(lambdaIdx);
end

save('lambda.mat')

