%%%% effect of lambda on procrastination 
rng shuffle
% increase lambda from lambda<1 to lambda>1, and then for lambda>1, increase lambda as well
NMonte = 1000;
deltas = 0.01;
k = 1;
Tmax = 10;
Tmin = 2;
TAll = randi([Tmin Tmax],1,NMonte);
alphaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
%JAll = zeros(1,NMonte);
JAll = 0.01*alphaAll.*rand(1,NMonte);
logk = -9+(1-(-9))*rand(1,NMonte);
k_DR = exp(logk);

Npara = 3;
lambdaAll = nan(NMonte,Npara);
lambdaAll(:,1) = rand(NMonte,1);
for nMonte = 1:NMonte
    loglambda = rand(1,Npara-1);
    lambdaVecend = sort(10.^(loglambda)); % lambda>1
    lambdaAll(nMonte,2:Npara) = lambdaVecend;
end

Ndaysofprocras = nan(NMonte,Npara);
Utask = nan(NMonte, Npara);
cost = nan(NMonte,Npara);
Nnetutility = nan(NMonte, Npara);
finalprop = nan(NMonte,Npara);

for nMonte=1:NMonte
    nMonte
        for idx = 1:Npara
        [temptemp,] = HyperbolicDiscountingFun(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte,idx),alphaAll(1,nMonte),betaAll(nMonte),k_DR(nMonte)); 
        
        if ~isempty(find(temptemp,1))
            Ndaysofprocras(nMonte,idx) = find(temptemp,1)-1;
        else
            Ndaysofprocras(nMonte,idx) = TAll(nMonte);
        end

        Utask(nMonte,idx)=alphaAll(nMonte)*(nansum(temptemp))^betaAll(nMonte);
        cost(nMonte,idx)=nansum(c1All(nMonte)*temptemp.^lambdaAll(nMonte,idx));
        Nnetutility(nMonte,idx) = Utask(nMonte,idx)-cost(nMonte,idx);  
        finalprop(nMonte,idx) = nansum(temptemp);

        end
            
end

%save('Hyperbolic_lambdaMonteCarlo_J0.mat')
%save('Hyperbolic_lambdaMonteCarlo_Jnonzero.mat')
deltadayspro = Ndaysofprocras(:,1)-Ndaysofprocras(:,2);
% find(deltadayspro<0) empty. 
deltafinalprop = finalprop(:,1)-finalprop(:,2);
% bidirectional
deltatotalcost = cost(:,1)-cost(:,2);
% bidirectional
b = Nnetutility(:,1)-Nnetutility(:,2);
% b bidirectional, mostly <=0, with a few exceptions.

%%%%%% when J>=0, conclusion when lambda changes from lambda<1 to lambda>1, days of
%%%%%% procrastination reduces, net utility increases or decreases.
%%%%%% task utility can increase or decrease or reserve, total cost can decrease
%%%%%% or increase or reserve.


deltadayspro = Ndaysofprocras(:,2)-Ndaysofprocras(:,3);
% J>0, bidirectional; J=0, deltadayspro>0
deltafinalprop = finalprop(:,2)-finalprop(:,3);
% J>0, bidirectional; J=0, deltafinalprop<0
deltatotalcost = cost(:,2)-cost(:,3);
% J>0, bidirectional; J=0, bidirectional
b = Nnetutility(:,2)-Nnetutility(:,3);
% J>0, bidirectional; J=0, b bidirectional, mostly <=0 with one exception

%%%%%% when J=0, conclusion: when lambda>1, when lambda increases, days of
%%%%%% procrastination reduces, net utility increases. task
%%%%%% utility can decrease or increase or reserve. cost can decrease or increase.

%%%%%% when J>0, conclusion in J=0 does not hold, days of procrastination
%%%%%% can be anything, net utility can be anything.

