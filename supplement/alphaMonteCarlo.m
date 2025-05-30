%%%%alpha increases, only in the delayed reward condition
NMonte = 10000;
k = 1; 
deltas = 0.01;

alphaAll = nan(2,NMonte);
alphaAll(1,:) = rand(1,NMonte);
alphaAll(2,:) = alphaAll(1,:)+(1-alphaAll(1,:)).*rand(1,NMonte);
Tmax=10;
TAll = randi([2 Tmax],1,NMonte);
gammaAll = rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll); % lambda>1
lambdaAll = 10*rand(1,NMonte);
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog;
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll(1,NMonte).*rand(1,NMonte);
%JAll = zeros(1,NMonte);

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
    [temp1,] = OptActStateSeq( deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(1,nMonte),betaAll(nMonte)); 
    WithoptActSeqMatrix1(nMonte,1:TAll(nMonte)) =  temp1;
    
    if ~isempty(find(temp1,1))
        Ndaysofprocras(nMonte,1) = find(temp1,1)-1;
    else
        Ndaysofprocras(nMonte,1) = TAll(nMonte);
    end
    [temp2,] = OptActStateSeq( deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(2,nMonte),betaAll(nMonte)); 
    WithoptActSeqMatrix2(nMonte,1:TAll(nMonte)) =  temp2;
    
    if ~isempty(find(temp2,1))
        Ndaysofprocras(nMonte,2) = find(temp2,1)-1;
    else
        Ndaysofprocras(nMonte,2) = TAll(nMonte);
    end
    
    U_task_1=alphaAll(1,nMonte)*(nansum(temp1))^betaAll(nMonte);
    cost_1=nansum(c1All(nMonte)*temp1.^lambdaAll(nMonte));
    Utotal_1 = U_task_1-cost_1;
    Nnetutility(nMonte,1) = Utotal_1;
    finalprop(nMonte,1) = nansum(temp1);
    U_task(nMonte,1) = U_task_1;
    cost(nMonte,1) = cost_1;
    
    U_task_2=alphaAll(2,nMonte)*(nansum(temp2))^betaAll(nMonte);
    cost_2=nansum(c1All(nMonte)*temp2.^lambdaAll(nMonte));
    Utotal_2 = U_task_2-cost_2;
    Nnetutility(nMonte,2) = Utotal_2; 
    finalprop(nMonte,2) = nansum(temp2);
    U_task(nMonte,2) = U_task_2;
    cost(nMonte,2) = cost_2;
    
end

%save('alphaMonteCarlo.mat')
deltadayspro = Ndaysofprocras(:,1)-Ndaysofprocras(:,2);
% % find(deltadayspro<0) % if no, days of procrastination is longer for smaller alpha
% % % Yes. 
deltafinalprop = finalprop(:,1)-finalprop(:,2);
%%%%%% find(deltafinalprop>0). if no, final proportion completed increases
%%%%%% over increasing alpha. Yes. 
deltatotalcost = cost(:,1)-cost(:,2);
%%%%%% find(deltatotalcost>0). If no, total cost increases over increasing
%%%%%% alpha. 
b = Nnetutility(:,1)-Nnetutility(:,2);
% % find(b>0) % if no, then net utility is larger for larger alpha
% % % Yes.always b<=0. 