%% c_1 increases, only in the delayed reward condition

NMonte  = 10000;
k = 1; 
deltas = 0.01;

c1All = nan(2,NMonte);
c1All(1,:) = 2*rand(1,NMonte);
c1All(2,:) = c1All(1,:)+ (2-c1All(1,:)).*rand(1,NMonte);

alphaAll = rand(1,NMonte);
Tmax = 10;
Tmin = 2;
TAll = randi([Tmin Tmax],1,NMonte);
gammaAll = rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll); % lambda>1
lambdaAll = 10*rand(1,NMonte);
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
JAll = 0.01*alphaAll.*rand(1,NMonte);
%JAll = zeros(1,NMonte);

Npara= 2;% 
Ndaysofprocras = nan(NMonte,Npara);
U_task = nan(NMonte,Npara);
cost = nan(NMonte,Npara);
Nnetutility = nan(NMonte, Npara);
finalprop = nan(NMonte,Npara);
WithoptActSeqMatrix1 = nan(NMonte,Tmax);
WithoptActSeqMatrix2 = nan(NMonte,Tmax);

for nMonte = 1: NMonte
    nMonte
    [temp1,] = OptActStateSeq( deltas,TAll(nMonte),k,JAll(nMonte),c1All(1,nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
    WithoptActSeqMatrix1(nMonte,1:TAll(nMonte)) =  temp1;
    if ~isempty(find(temp1,1))
        Ndaysofprocras(nMonte,1) = find(temp1,1)-1;
    else
        Ndaysofprocras(nMonte,1) = TAll(nMonte);
    end
    [temp2,] = OptActStateSeq( deltas,TAll(nMonte),k,JAll(nMonte),c1All(2,nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
    WithoptActSeqMatrix2(nMonte,1:TAll(nMonte)) =  temp2;
    if ~isempty(find(temp2,1))
        Ndaysofprocras(nMonte,2) = find(temp2,1)-1;
    else
        Ndaysofprocras(nMonte,2) = TAll(nMonte);
    end
    
    U_task(nMonte,1)=alphaAll(nMonte)*(nansum(temp1))^betaAll(nMonte);
    cost(nMonte,1)=nansum(c1All(1,nMonte)*temp1.^lambdaAll(nMonte));
    Nnetutility(nMonte,1) = U_task(nMonte,1)-cost(nMonte,1);
    finalprop(nMonte,1) = nansum(temp1);
    
    U_task(nMonte,2)=alphaAll(nMonte)*(nansum(temp2))^betaAll(nMonte);
    cost(nMonte,2)=nansum(c1All(2,nMonte)*temp2.^lambdaAll(nMonte));
    Nnetutility(nMonte,2) = U_task(nMonte,2)-cost(nMonte,2);
    finalprop(nMonte,2) = nansum(temp2);
    
end

%save('CmaxMonteCarlo_J0.mat') % this is for the case when J=0
%save('CmaxMonteCarlo_Jnonzero.mat') % this is for the case when J>0
deltadayspro = Ndaysofprocras(:,1)-Ndaysofprocras(:,2);
% find(deltadayspro>0) if no, days of procrastination is longer for higher C_max 
% Yes. 
deltafinalprop = finalprop(:,1)-finalprop(:,2);
%%%%%% find(deltafinalprop<0). If no, final proportion completed decreases
%%%%%% for higher C_max. 
deltatotalcost = cost(:,1)-cost(:,2);
%%%%%% find(deltatotalcost>0). total cost is bidirectional
b = Nnetutility(:,1)-Nnetutility(:,2);
% find(b<0) 

%%%%%%% conclusion: when J=0, increasing c_max results in: procrastinate
%%%%%%% for longer,  final proportion completed decrease, and cost can
%%%%%%% increase or decrease, net utility decrease.
%%%%%%% when J>0, change in everything is bidirectional (days of procrastination, net utility, task utility or cost)
