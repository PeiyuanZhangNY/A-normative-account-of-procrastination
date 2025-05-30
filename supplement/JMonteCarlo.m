%%%%J increases, only in the delayed reward condition

NMonte  = 10000;
k = 1; 
deltas = 0.01;

Tmax=10;
TAll = randi([2 Tmax],1,NMonte);
alphaAll = rand(1,NMonte);
gammaAll = rand(1,NMonte);
lambdaAll = 10*rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll); % lambda>1
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog;
c1All = 2*rand(1,NMonte);
JAll = nan(2,NMonte);
JAll(1,:) = 0.01*alphaAll.*rand(1,NMonte);
JAll(2,:) = JAll(1,:)+(0.01*alphaAll-JAll(1,:))*rand;

Npara= 2;% 
Ndaysofprocras = nan(NMonte,Npara);
Nnetutility = nan(NMonte, Npara);
finalprop = nan(NMonte,Npara);
cost = nan(NMonte,Npara);
WithoptActSeqMatrix1 = nan(NMonte,Tmax);
WithoptActSeqMatrix2 = nan(NMonte,Tmax);

for nMonte = 1: NMonte
    nMonte
    [temp1,] = OptActStateSeq( deltas,TAll(nMonte),k,JAll(1,nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
    WithoptActSeqMatrix1(nMonte,1:TAll(nMonte)) =  temp1;
    if ~isempty(find(temp1,1))
        Ndaysofprocras(nMonte,1) = find(temp1,1)-1;
    else
        Ndaysofprocras(nMonte,1) = TAll(nMonte);
    end
    
    [temp2,] = OptActStateSeq( deltas,TAll(nMonte),k,JAll(2,nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
    WithoptActSeqMatrix2(nMonte,1:TAll(nMonte)) =  temp2;
    if ~isempty(find(temp2,1))
        Ndaysofprocras(nMonte,2) = find(temp2,1)-1;
    else
        Ndaysofprocras(nMonte,2) = TAll(nMonte);
    end
    
    U_task_1=alphaAll(nMonte)*(nansum(temp1))^betaAll(nMonte);
    cost_1=nansum(c1All(nMonte)*temp1.^lambdaAll(nMonte));
    Utotal_1 = U_task_1-cost_1;
    finalprop(nMonte,1) = nansum(temp1);
    Nnetutility(nMonte,1) = Utotal_1;
    cost(nMonte,1) = cost_1;
    
    U_task_2=alphaAll(nMonte)*(nansum(temp2))^betaAll(nMonte);
    cost_2=nansum(c1All(nMonte)*temp2.^lambdaAll(nMonte));
    Utotal_2 = U_task_2-cost_2;
    Nnetutility(nMonte,2) = Utotal_2;   
    finalprop(nMonte,2) = nansum(temp2);
    cost(nMonte,2) = cost_2;
end

%save('JMonteCarlo.mat')

deltadayspro = Ndaysofprocras(:,1)-Ndaysofprocras(:,2);
%find(deltadayspro<0) 
% % Yes. 
deltafinalprop = finalprop(:,1)-finalprop(:,2);
%find(deltafinalprop>0)
% % yes. 
deltatotalcost = cost(:,1)-cost(:,2);
%find(deltatotalcost>0) bidirectional 
b = Nnetutility(:,1)-Nnetutility(:,2);
find(b<0) 
% % Yes.b is always>=0, meaning net utility decreases when J increases.
