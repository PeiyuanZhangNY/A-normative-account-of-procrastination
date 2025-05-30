NMonte  = 10000;
k = 1; 
deltas = 0.01;
TAll = nan(2,NMonte);
Tmax=10;
TAll(1,:) = randi([2 Tmax-1],1,NMonte);
TAll(2,:) = TAll(1,:)+1; % avoid when T does not increase by one unit, it is possible T1 = 2, [0.47,0.53], T2=10, [0,0,0,0,0,0,0.21,0.24,0.26,0.29], actually between T1 and T2, there is case [0.24,0.26,0.29]
% for nMonte = 1:NMonte
%     TAll(2,nMonte) = randi([TAll(1,nMonte)+1 Tmax]);
% end
gammaAll = rand(1,NMonte);
loglambdaAll = rand(1,NMonte); 
lambdaAll = 10.^(loglambdaAll); % lambda>1
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog;  % 0.1<beta<10 so beta can be <lambda, or >lambda
alphaAll = rand(1,NMonte);
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll.*rand(1,NMonte);
Npara= 2;% 
Ndaysofprocras = nan(NMonte,Npara);
Nnetutility = nan(NMonte, Npara);
finalprop = nan(NMonte,Npara);
U_task = nan(NMonte,Npara);
cost = nan(NMonte,Npara);
WithoptActSeqMatrix1 = nan(NMonte,Tmax);
WithoptActSeqMatrix2 = nan(NMonte,Tmax);

allzerosmark1 = nan(NMonte,Tmax); % cases of [0,0]
allzerosmark2 = nan(NMonte,Tmax); % cases of [0,0,0]
repeatZeroBegin1 = nan(NMonte,Tmax); % cases of [0,0.2,0.3]
repeatZeroBegin2 = nan(NMonte,Tmax); % when T increases by 1, then [0,0,0.2,0.3], the sequence is repeated

for nMonte = 1: NMonte
    nMonte
    [temp1,] = OptActStateSeq( deltas,TAll(1,nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
    [temp2,] = OptActStateSeq( deltas,TAll(2,nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
    if ~isempty(find(temp1,1))
        temp1without0 = temp1(find(temp1,1):end);
    end
    
    if ~isempty(find(temp2,1))
        temp2without0 = temp2(find(temp2,1):end);
    end
    
    if sum(temp1)+sum(temp2)>0 % we don't count those days temp1=[0,0], temp2=[0,0,0], when the normalized procrastination days is naively 1 and unchanged with increased T.
        if ~isequal(temp1without0,temp2without0)     % delete the cases where the policy get procrastination because of
        % deltas resolution. for example T=3, [0.01,0.02,0.05]. T=4, [0,0.01,0.02,0.05]
        % time is too long then necessary. Because in our initial assumption, effort a is a
        % continuous variable, so it is possible a-->0 when T is infinite. but
        % here as the resolution of a is deltas. 
            WithoptActSeqMatrix1(nMonte,1:TAll(1,nMonte)) = temp1;
            if ~isempty(find(temp1,1))
                Ndaysofprocras(nMonte,1) = find(temp1,1)-1;
            else
                Ndaysofprocras(nMonte,1) = TAll(1,nMonte);
            end

            WithoptActSeqMatrix2(nMonte,1:TAll(2,nMonte)) = temp2;
            if ~isempty(find(temp2,1))
                Ndaysofprocras(nMonte,2) = find(temp2,1)-1;
            else
                Ndaysofprocras(nMonte,2) = TAll(2,nMonte);
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
        else
            repeatZeroBegin1(nMonte,1:TAll(1,nMonte)) = temp1;
            repeatZeroBegin2(nMonte,1:TAll(2,nMonte)) = temp2;
        end
    else
        allzerosmark1(nMonte,1:TAll(1,nMonte)) = temp1;
        allzerosmark2(nMonte,1:TAll(2,nMonte)) = temp2;
    end
end

a = Ndaysofprocras(:,1)-Ndaysofprocras(:,2);
% find(a<0) % if no, then days of procrastination is longer
% % Yes. 
b = Nnetutility(:,1)-Nnetutility(:,2);
% find(b>0) % if no, then net utility in the lambda>1 condition larger than the delayed condition
% % Yes.b always<=0.
normalizedNdaysofprocras = nan(NMonte,2);
normalizedNdaysofprocras(:,1) = Ndaysofprocras(:,1)./(TAll(1,:)');
normalizedNdaysofprocras(:,2) = Ndaysofprocras(:,2)./(TAll(2,:)');

figure
for nMonte =1:NMonte
    plot(1:2,normalizedNdaysofprocras(nMonte,:),'ko-')
    hold on
end
% normalized days of procrastination switches from 1 to 0 with increased T.
% 
figure
for nMonte =1:NMonte
    plot(1:2,finalprop(nMonte,:),'ko-')
    hold on
end
% finalprop increases with T.
figure
for nMonte =1:NMonte
    plot(1:2,cost(nMonte,:),'ko-')
    hold on
end
% bidirectional 
figure
for nMonte =1:NMonte
    plot(1:2,Nnetutility(nMonte,:),'ko-')
    hold on
end
% save('TMonteCarlo.mat')
%%%%%%%% conclusion: when we increase T, first from [0,0] to [0,0,0], then to work for all days [0.3,0.4,0.4], 
%%%%%%% then still work for all days [0.1, 0.1, 0.2,0.3,0.3]
%%%%%%% if we calculate the normalized days of procrastination, change from
%%%%%%% 1 to 0.
