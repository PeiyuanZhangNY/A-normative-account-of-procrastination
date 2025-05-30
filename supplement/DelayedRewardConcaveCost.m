close all;
% Monte Carlo random sampling 
NMonte  = 10000;
deltas = 0.01;
k = 1;
Tmax = 10;
Tmin = 2;
TAll = randi([Tmin Tmax],1,NMonte);
alphaAll = rand(1,NMonte); 
gammaAll = rand(1,NMonte);
c1All = 2*rand(1,NMonte);
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog;
JAll = 0.01*alphaAll.*rand(1,NMonte);
lambdaAll = rand(1,NMonte); 

WithoptActSeqMatrix=nan(NMonte,Tmax);
for nMonte = 1:NMonte
    nMonte
    [WithoptActSeqMatrix(nMonte,1:TAll(nMonte)), ]=OptActStateSeq( deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte));    
end

%% check whether last minte worker conclusion holds true
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

