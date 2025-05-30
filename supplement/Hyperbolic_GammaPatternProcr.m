%%%%%%%%%% hyperbolic discounting, simulations for the agent works little in the early days, and work more and more towards the deadline. 
NMonte  = 10000;
k = 1;
Tmax = 10;
Tmin = 2;
TAll = randi([Tmin Tmax],1,NMonte);
alphaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
loglambdaAll = rand(1,NMonte); 
lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
%JAll = zeros(1,NMonte);
JAll = 0.01*alphaAll.*rand(1,NMonte);

logk = -9+(1-(-9))*rand(1,NMonte);
k_DR = exp(logk);

WithoptActSeqMatrix=nan(NMonte,Tmax);
for nMonte = 1:NMonte
    nMonte
    [WithoptActSeqMatrix(nMonte,1:TAll(nMonte)), ]=HyperbolicDiscountingFun(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DR(nMonte)); 
end

% draw figure to show the raw data 
figure
for nMonte = 1:100%NMonte
    plot(1:Tmax,WithoptActSeqMatrix(nMonte,:),'bo-')
    hold on
end 
xlabel('time t')
xlim([1,Tmax])
ylabel('effort a')

%% test whether the conclusion (work more and more towards the deadline) hold true
tplus1minust = nan(1,NMonte*Tmax);
n=1; 
for nMonte = 1:NMonte
     effortvec = WithoptActSeqMatrix(nMonte,1:TAll(nMonte));
     for effortidx = 1:length(effortvec)-1
         tplus1minust(1,n)=effortvec(effortidx+1)-effortvec(effortidx);
         if tplus1minust(1,n)<0
             disp(nMonte)
         end
             
         n=n+1;
     end
end
tplus1minust(tplus1minust<0)
figure
histogram(tplus1minust)

