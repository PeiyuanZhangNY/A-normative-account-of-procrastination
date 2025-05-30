%%%%%%%% hyperbolic discounting, simulations for parameter space for the pattern ?working at the last minute?
NMonte  = 10000;
k = 1;
Tmax = 10;
Tmin = 2;
TAll = randi([Tmin Tmax],1,NMonte);
alphaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog;
JAll = 0.01*alphaAll.*rand(1,NMonte);
lambdaAll = rand(1,NMonte); 

logk = -9+(1-(-9))*rand(1,NMonte);
k_DR = exp(logk);

WithoptActSeqMatrix=nan(NMonte,Tmax);
for nMonte = 1:NMonte
    nMonte
    [WithoptActSeqMatrix(nMonte,1:TAll(nMonte)), ]=HyperbolicDiscountingFun(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DR(nMonte));    
end

%% check whether last minute worker conclusion holds true
numberofnonezeros = nan(1,NMonte);
for nMonte=1:NMonte
    numberofnonezeros(nMonte)=nnz(WithoptActSeqMatrix(nMonte,1:TAll(nMonte)));
    indexnonzero = find(WithoptActSeqMatrix(nMonte,1:TAll(nMonte))>0);
    if ~isempty(indexnonzero)
        if indexnonzero~=TAll(nMonte)
            fprintf(num2str(nMonte))
        end
    end
end
unique(numberofnonezeros) %0, 1 % yes last minute worker conclusion holds true 

