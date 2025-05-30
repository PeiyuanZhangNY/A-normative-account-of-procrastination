%%%%%%%%%%% the agent works little in the early days, and work more and more towards the deadline.
NMonte  = 10000;
deltas = 0.01;
k = 1;
Tmax = 10;
Tmin = 2;
TAll = randi([Tmin Tmax],1,NMonte);
alphaAll = rand(1,NMonte); 
gammaAll = rand(1,NMonte); 
c1All = 2*rand(1,NMonte);
loglambdaAll = rand(1,NMonte); 
lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
%JAll = zeros(1,NMonte);
JAll = 0.01*alphaAll.*rand(1,NMonte);

WithoptActSeqMatrix=nan(NMonte,Tmax);
for nMonte = 1:NMonte
    %nMonte
    [WithoptActSeqMatrix(nMonte,1:TAll(nMonte)), ]=OptActStateSeq(deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
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
    nMonte
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

% the effort on the last day in the special case gamma=0
% lasteffort = nan(NMonte,1);
% for nMonte = 1:NMonte
%     lasteffort(nMonte,1) = WithoptActSeqMatrix(nMonte,TAll(nMonte));
% end

% %% test whether the conclusion (\gamma=1, work steadily) holds true
% %%% comment off gammaAll = rand(1,NMonte); and write instead gammaAll = ones(1,NMonte);
% tplus1minust = nan(1,NMonte*Tmax);
% n=1; 
% for nMonte = 1:NMonte
%      effortvec = WithoptActSeqMatrix(nMonte,1:TAll(nMonte));
%      for effortidx = 1:length(effortvec)-1
%          tplus1minust(1,n)=effortvec(effortidx+1)-effortvec(effortidx);
%          if tplus1minust(1,n)>0 || tplus1minust(1,n)<0
%              disp(nMonte)
%          end
%              
%          n=n+1;
%      end
% end
% 
% figure
% histogram(tplus1minust)
% 
% figure
% histogram([tplus1minust(tplus1minust<0),tplus1minust(tplus1minust>0)])

