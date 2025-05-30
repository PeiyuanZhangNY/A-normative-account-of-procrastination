% test whether the number of days of procratination decreases over increasing gamma
% original file DelayedRewardGammaDaysUtility.m 
NMonte  = 10000;
deltas = 0.01;
k = 1;
Tmax = 10;
Tmin = 2;
TAll = randi([Tmin Tmax],1,NMonte);
alphaAll = rand(1,NMonte);  
c1All = 2*rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
lambdaAll = 10*rand(1,NMonte);
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
%JAll = zeros(1,NMonte);
JAll = 0.01*alphaAll.*rand(1,NMonte);
k_DRAll = nan(2,NMonte); 
logk = -9+(1-(-9))*rand(1,NMonte);
k_DR = exp(logk);
k_DRAll(1,:) = k_DR;
k_DRAll(2,:) = k_DRAll(1,:) + (exp(1)-k_DRAll(1,:)).*rand(1,NMonte);

Npara= 2;% 
Ndaysofprocras = nan(NMonte,Npara);
Nnetutility = nan(NMonte, Npara);
finalprop = nan(NMonte,Npara);
U_task = nan(NMonte,Npara);
cost = nan(NMonte,Npara);
WithoptActSeqMatrix1 = nan(NMonte,Tmax);
WithoptActSeqMatrix2 = nan(NMonte,Tmax);

for nMonte = 1: NMonte
    nMonte
    [temp1,] = HyperbolicDiscountingFun(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DRAll(1,nMonte)); 
    WithoptActSeqMatrix1(nMonte,1:TAll(nMonte)) =  temp1;
    
    if ~isempty(find(temp1,1))
        Ndaysofprocras(nMonte,1) = find(temp1,1)-1;
    else
        Ndaysofprocras(nMonte,1) = TAll(nMonte);
    end
    [temp2,] = HyperbolicDiscountingFun(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DRAll(2,nMonte)); 
    WithoptActSeqMatrix2(nMonte,1:TAll(nMonte)) =  temp2;
    
    if ~isempty(find(temp2,1))
        Ndaysofprocras(nMonte,2) = find(temp2,1)-1;
    else
        Ndaysofprocras(nMonte,2) = TAll(nMonte);
    end
    
    U_task_1=alphaAll(nMonte)*(nansum(temp1))^betaAll(nMonte);
    cost_1=nansum(c1All(nMonte)*temp1.^lambdaAll(nMonte));
    Utotal_1 = U_task_1-cost_1;
    Nnetutility(nMonte,1) = Utotal_1;
    finalprop(nMonte,1) = nansum(temp1);
    U_task(nMonte,1) = U_task_1;
    cost(nMonte,1) = cost_1;
    
    U_task_2=alphaAll(nMonte)*(nansum(temp2))^betaAll(nMonte);
    cost_2=nansum(c1All(nMonte)*temp2.^lambdaAll(nMonte));
    Utotal_2 = U_task_2-cost_2;
    Nnetutility(nMonte,2) = Utotal_2; 
    finalprop(nMonte,2) = nansum(temp2);
    U_task(nMonte,2) = U_task_2;
    cost(nMonte,2) = cost_2;
    
end

%save('Hyperbolic_k_DRMonteCarlo.mat')
deltadayspro = Ndaysofprocras(:,1)-Ndaysofprocras(:,2);
% % find(deltadayspro>0) % if no, days of procrastination is longer for
% larger k
% % % Yes. 
deltafinalprop = finalprop(:,1)-finalprop(:,2);
%%%%%% find(deltafinalprop<0). if no, final proportion completed increases
%%%%%% over increasing k. Yes. 
deltatotalcost = cost(:,1)-cost(:,2);
%%%%%% find(deltatotalcost<0). total cost is bidirectional.
b = Nnetutility(:,1)-Nnetutility(:,2);
%%%%%% b is mostly >=0, except 1 case. 

