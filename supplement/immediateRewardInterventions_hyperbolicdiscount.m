close all;
clear 
rng('shuffle')
% data saved in 'ImmediateRewardIntervention_2.mat'
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);
tic

deltas=0.01;
k=1;
NMonte = 10000;

alphaAll = rand(1,NMonte);
%Tmax=14;
%TAll = randi([6,Tmax],1,NMonte);
TAll = 10*ones(1,NMonte); %fixed at 10 days. 
logk = -9+(1-(-9))*rand(1,NMonte);
k_DR = exp(logk);
loglambdaAll = rand(1,NMonte); 
lambdaAll = 10.^(loglambdaAll); % lambda>1
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog;
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll(1,NMonte).*rand(1,NMonte);
N_milestone = randi([2,5],1,NMonte); % has to be smaller than TAll

Npara= 4;% delayed reward and 3 immediate reward interventions 
Ndaysofprocras = nan(NMonte,Npara);
MUCD = nan(NMonte,Npara); % mean unit completion day
%Normalizeddaysofprocras = nan(NMonte,Npara);
%Nnetutility = nan(NMonte, Npara);
finisheday = nan(NMonte,Npara);
finalprop = nan(NMonte,Npara);

WithoptActSeqMatrix1 = nan(NMonte,TAll(1)); % delayed reward
WithoptActSeqMatrix2 = nan(NMonte,TAll(1)); % reward upon task completion 
WithoptActSeqMatrix3 = nan(NMonte,TAll(1)); % reward each milestone
WithoptActSeqMatrix4 = nan(NMonte,TAll(1)); % reward each unit of progress

for nMonte = 1588: NMonte
    nMonte
    
    [temp1,] = HyperbolicDiscountingFun(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DR(nMonte)); 
    WithoptActSeqMatrix1(nMonte,1:TAll(nMonte)) =  temp1;
    
    if ~isempty(find(temp1,1))
        Ndaysofprocras(nMonte,1) = find(temp1,1)-1;
    else
        Ndaysofprocras(nMonte,1) = TAll(nMonte);
    end
    finalprop(nMonte,1) = nansum(temp1);
    if nansum(temp1)==1
        MUCD(nMonte,1) = nansum(temp1'.*(1:TAll(nMonte)));
        finisheday(nMonte,1) = find(cumsum(temp1)==1,1);
    end
    
    
    
    [temp2,] = HyperbolicDiscountingFun_RewardCompletionTime(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DR(nMonte)); 
    WithoptActSeqMatrix2(nMonte,1:TAll(nMonte)) =  temp2;
    
    if ~isempty(find(temp2,1))
        Ndaysofprocras(nMonte,2) = find(temp2,1)-1;
    else
        Ndaysofprocras(nMonte,2) = TAll(nMonte);
    end
  
    finalprop(nMonte,2) = nansum(temp2);
    
    if nansum(temp2)==1
        MUCD(nMonte,2) = nansum(temp2'.*(1:TAll(nMonte)));
        finisheday(nMonte,2) = find(cumsum(temp2)==1,1);
    end
    
    [temp3,] = HyperbolicDiscountingFun_RewardMilestons(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DR(nMonte),N_milestone(nMonte)); 
    WithoptActSeqMatrix3(nMonte,1:TAll(nMonte)) =  temp3;
    
    if ~isempty(find(temp3,1))
        Ndaysofprocras(nMonte,3) = find(temp3,1)-1;
    else
        Ndaysofprocras(nMonte,3) = TAll(nMonte);
    end
  
    finalprop(nMonte,3) = nansum(temp3);
    
    if nansum(temp3)==1
        MUCD(nMonte,3) = nansum(temp3'.*(1:TAll(nMonte)));
        finisheday(nMonte,3) = find(cumsum(temp3)==1,1);
    end    

    [temp4,] = HyperbolicDiscountingFun_RewardEachUnitOfProgress(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DR(nMonte)); 
    WithoptActSeqMatrix4(nMonte,1:TAll(nMonte)) =  temp4;
    
    if ~isempty(find(temp4,1))
        Ndaysofprocras(nMonte,4) = find(temp4,1)-1;
    else
        Ndaysofprocras(nMonte,4) = TAll(nMonte);
    end
  
    finalprop(nMonte,4) = nansum(temp4);
    
    if nansum(temp4)==1
        MUCD(nMonte,4) = nansum(temp4'.*(1:TAll(nMonte)));
        finisheday(nMonte,4) = find(cumsum(temp4)==1,1);
    end 
    
end
toc % 2 minutes for NMonte=1000;  15 minutes for NMonte=10000;

% calculate all utilities
U_task_with = nan(NMonte,Npara);
cost_with = nan(NMonte,Npara);
Utotal_with = nan(NMonte,Npara);

for nMonte = 1:NMonte
    U_task_with(nMonte,1) = alphaAll(nMonte)*(k*nansum(WithoptActSeqMatrix1(nMonte,:))^betaAll(nMonte));
    U_task_with(nMonte,2) = alphaAll(nMonte)*(k*nansum(WithoptActSeqMatrix2(nMonte,:))^betaAll(nMonte));
    U_task_with(nMonte,3) = alphaAll(nMonte)*(k*nansum(WithoptActSeqMatrix3(nMonte,:))^betaAll(nMonte));
    U_task_with(nMonte,4) = alphaAll(nMonte)*(k*nansum(WithoptActSeqMatrix4(nMonte,:))^betaAll(nMonte));
    
    cost_with(nMonte,1) = nansum(c1All(nMonte)*WithoptActSeqMatrix1(nMonte,:).^lambdaAll(nMonte));
    cost_with(nMonte,2) = nansum(c1All(nMonte)*WithoptActSeqMatrix2(nMonte,:).^lambdaAll(nMonte));
    cost_with(nMonte,3) = nansum(c1All(nMonte)*WithoptActSeqMatrix3(nMonte,:).^lambdaAll(nMonte));
    cost_with(nMonte,4) = nansum(c1All(nMonte)*WithoptActSeqMatrix4(nMonte,:).^lambdaAll(nMonte));
    
    Utotal_with(nMonte,1) = U_task_with(nMonte,1)-cost_with(nMonte,1);
    Utotal_with(nMonte,2) = U_task_with(nMonte,2)-cost_with(nMonte,2);
    Utotal_with(nMonte,3) = U_task_with(nMonte,3)-cost_with(nMonte,3);
    Utotal_with(nMonte,4) = U_task_with(nMonte,4)-cost_with(nMonte,4);
end

save('ImmediateRewardInterventions_hyperbolicdiscount.mat')

%% check if the average result hold true for each simulation
% days of prorastination delayed>upon task completion>upon milestones>upon
% unit of work
conditionMet = Ndaysofprocras(:,1) >= Ndaysofprocras(:,2) & ...
               Ndaysofprocras(:,2) >= Ndaysofprocras(:,3) & ...
               Ndaysofprocras(:,3) >= Ndaysofprocras(:,4);

% Find the indices of rows where the condition is violated
rowsViolating = find(~conditionMet);
% result: there are a few violations

% check if immediate reward inteventions reduces procrastination, (not necessarily with higher level of immediacy, fewer days of procrastination)
conditionMet = Ndaysofprocras(:,1) >= Ndaysofprocras(:,2) & ...
               Ndaysofprocras(:,1) >= Ndaysofprocras(:,3) & ...
               Ndaysofprocras(:,1) >= Ndaysofprocras(:,4);
rowsViolating = find(~conditionMet);
% result: YES! 

% task completion day 
rowsWithoutNaN = all(~isnan(finisheday), 2);
conditionMet = (finisheday(rowsWithoutNaN,1) >= finisheday(rowsWithoutNaN,2)) & ...
               (finisheday(rowsWithoutNaN,1) >= finisheday(rowsWithoutNaN,3)) & ...
               (finisheday(rowsWithoutNaN,1) >= finisheday(rowsWithoutNaN,4));
% Get the original row indices (before removing NaNs)
originalRowIndices = find(rowsWithoutNaN);

% Find the indices of the original rows that violate the condition
rowsViolatingCondition = originalRowIndices(~conditionMet);
% result: 0×1 empty; YES! In terms of task completion day, in contrast to the delayed reward condition, where agents always complete tasks on the last day, in every single simulation, offering immediate rewards helps agents complete the task ahead of the deadline. There is no monotonic relationship between task completion day and the immediacy level of interventions.

% final proporiton completed
conditionMet = finalprop(:,1) <= finalprop(:,2) & ...
               finalprop(:,2) <= finalprop(:,3) & ...
               finalprop(:,3) <= finalprop(:,4);

% Find the indices of rows where the condition is violated
rowsViolating = find(~conditionMet);
% result: NO! In terms of the final proportion completed, on average, offering immediate reward helps agents get more work done. The higher the immediacy level of interventions, the higher the final proportion completed. However, it does not always hold for every single simulation.

conditionMet = finalprop(:,1) <= finalprop(:,2) & ...
               finalprop(:,1) <= finalprop(:,3) & ...
               finalprop(:,1) <= finalprop(:,4);

% Find the indices of rows where the condition is violated
rowsViolating = find(~conditionMet);
% result: NO!


% final reward
conditionMet = U_task_with(:,1) <= U_task_with(:,2) & ...
               U_task_with(:,2) <= U_task_with(:,3) & ...
               U_task_with(:,3) <= U_task_with(:,4);

rowsViolating = find(~conditionMet);
% result: NO!

% total cost
conditionMet = cost_with(:,1) <= cost_with(:,2) & ...
               cost_with(:,1) <= cost_with(:,3) & ...
               cost_with(:,1) <= cost_with(:,4);

rowsViolating = find(~conditionMet);
% result: NO!

% net utility
conditionMet = Utotal_with(:,1) <= Utotal_with(:,2) & ...
               Utotal_with(:,1) <= Utotal_with(:,3) & ...
               Utotal_with(:,1) <= Utotal_with(:,4);

rowsViolating = find(~conditionMet);

%% mean of all the outcomes

figure
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

subplot(1,4,1)
b=bar([1,2,3,4],nanmean(Ndaysofprocras,1),'FaceColor','flat');
b.CData(1,:) = [0,0,0];
b.CData(2,:) = greencolor(6,:);
b.CData(3,:) = greencolor(4,:);
b.CData(4,:) = greencolor(2,:);
b.EdgeColor = 'flat';
xticks([1,2,3,4])
%xticklabels({'delayed reward','reward upon task completion','reward each milestone (n\in [2,5])','reward each unit of proress'})
xticklabels({'delayed reward','reward upon task completion','reward each milestone','reward each unit of proress'})
yticks(0:2:TAll(1));
xtickangle(45);
ylabel('average days of procrastination')
ylim([0,TAll(1)])
box off
set(gca,'TickDir','out');

subplot(1,4,4)
% four groups of three bars
threecolor = [27,158,119;217,95,2;117,112,179]/255;
meanU_task = nanmean(U_task_with,1);
meanCost = nanmean(cost_with,1);
meanUnet = nanmean(Utotal_with,1);
data = [meanU_task(1),meanCost(1),meanUnet(1);meanU_task(2),meanCost(2),meanUnet(2);meanU_task(3),meanCost(3),meanUnet(3);meanU_task(4),meanCost(4),meanUnet(4)];
b = bar(data, 'grouped');

% Set the colors for each bar
for i = 1:size(data, 2)  % Loop over the number of bars in each group
    b(i).FaceColor = 'flat';
    b(i).CData = repmat(threecolor(i, :), size(data, 1), 1);
    b(i).EdgeColor = 'flat';
end
legend({'final reward', 'total cost', 'net utility'})
legend boxoff
xticks([1,2,3,4])
xticklabels({'delayed reward','reward upon task completion','reward each milestone','reward each unit of proress'})
yticks(0:0.2:1);
xtickangle(45);
ylabel('utility')
ylim([0,1])
box off
set(gca,'TickDir','out');

subplot(1,4,2)
b=bar([1,2,3,4],nanmean(finisheday,1),'FaceColor','flat');
b.CData(1,:) = [0,0,0];
b.CData(2,:) = greencolor(6,:);
b.CData(3,:) = greencolor(4,:);
b.CData(4,:) = greencolor(2,:);
b.EdgeColor = 'flat';
xticks([1,2,3,4])
xticklabels({'delayed reward','reward upon task completion','reward each milestone','reward each unit of proress'})
yticks(0:2:TAll(1));
xtickangle(45);
ylabel('average task completion day')
ylim([0,TAll(1)])
box off
set(gca,'TickDir','out');

subplot(1,4,3)
b=bar([1,2,3,4],nanmean(finalprop,1),'FaceColor','flat');
b.CData(1,:) = [0,0,0];
b.CData(2,:) = greencolor(6,:);
b.CData(3,:) = greencolor(4,:);
b.CData(4,:) = greencolor(2,:);
b.EdgeColor = 'flat';
xticks([1,2,3,4])
xticklabels({'delayed reward','reward upon task completion','reward each milestone','reward each unit of proress'})
yticks(0:0.2:1);
xtickangle(45);
ylabel('average final proportion completed')
ylim([0,1])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=300;
set(gcf,'position',[x0,y0,width,height])




