close all;
clear 
rng('shuffle')
% data saved in 'ImmediateRewardIntervention_2.mat'
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);
tic
T = 21;
N_milestone = 3;
deltas=0.01;
k=1;
NMonte = 10000;

alphaAll = rand(1,NMonte);
%Tmax=14;
%TAll = randi([6,Tmax],1,NMonte);
TAll = T*ones(1,NMonte);  
logk = -9+(1-(-9))*rand(1,NMonte);
k_DR = exp(logk);
loglambdaAll = rand(1,NMonte); 
lambdaAll = 10.^(loglambdaAll); % lambda>1
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog;
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll(1,NMonte).*rand(1,NMonte);


Npara= 2;% delayed reward and interim deadline
Ndaysofprocras = nan(NMonte,Npara);
MUCD = nan(NMonte,Npara); % mean unit completion day
%Normalizeddaysofprocras = nan(NMonte,Npara);
%Nnetutility = nan(NMonte, Npara);
finisheday = nan(NMonte,Npara);
finalprop = nan(NMonte,Npara);

WithoptActSeqMatrix1 = nan(NMonte,TAll(1)); % delayed reward
WithoptActSeqMatrix2 = nan(NMonte,TAll(1)); % interim deadline

net_earning_all = nan(NMonte,1);
%penalty_unit = 0.03;
penalty_unit_all = 0.1*rand(1,NMonte);

for nMonte = 1: NMonte
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
    
    
    
    [temp2,~,net_earning] = HyperbolicDiscountingFun_InterimDeadline(TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),k_DR(nMonte),N_milestone,penalty_unit_all(nMonte)); 
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
    net_earning_all(nMonte) = net_earning;
    
end
toc 

% calculate all utilities
U_task_with = nan(NMonte,Npara);
cost_with = nan(NMonte,Npara);
Utotal_with = nan(NMonte,Npara);

for nMonte = 1:NMonte
    U_task_with(nMonte,1) = alphaAll(nMonte)*(k*nansum(WithoptActSeqMatrix1(nMonte,:))^betaAll(nMonte));
    U_task_with(nMonte,2) = alphaAll(nMonte)*(k*nansum(WithoptActSeqMatrix2(nMonte,:))^betaAll(nMonte));
    
    cost_with(nMonte,1) = nansum(c1All(nMonte)*WithoptActSeqMatrix1(nMonte,:).^lambdaAll(nMonte));
    cost_with(nMonte,2) = nansum(c1All(nMonte)*WithoptActSeqMatrix2(nMonte,:).^lambdaAll(nMonte));
    
    Utotal_with(nMonte,1) = U_task_with(nMonte,1)-cost_with(nMonte,1);
    Utotal_with(nMonte,2) = U_task_with(nMonte,2)-cost_with(nMonte,2);
end

%save('interim_deadline_revised_hyperbolicdiscount.mat')

save('InterimDeadlineInterventions_hyperbolicdiscount.mat')
%% check if the average result hold true for each simulation
% days of prorastination delayed>upon task completion>upon milestones>upon
% unit of work
conditionMet = Ndaysofprocras(:,1) >= Ndaysofprocras(:,2);
% Find the indices of rows where the condition is violated
rowsViolating = find(~conditionMet);
% result: 0×1 empty; YES!  In terms of days of procrastination, in every single simulation, offering immediate rewards reduces the number of days of delay. The higher the immediacy level of interventions, the fewer days of procrastination.

% task completion day 
rowsWithoutNaN = all(~isnan(finisheday), 2);
conditionMet = finisheday(rowsWithoutNaN,1) >= finisheday(rowsWithoutNaN,2);
% Get the original row indices (before removing NaNs)
originalRowIndices = find(rowsWithoutNaN);

% Find the indices of the original rows that violate the condition
rowsViolatingCondition = originalRowIndices(~conditionMet);
% result: 0×1 empty; YES! In terms of task completion day, in contrast to the delayed reward condition, where agents always complete tasks on the last day, in every single simulation, offering immediate rewards helps agents complete the task ahead of the deadline. There is no monotonic relationship between task completion day and the immediacy level of interventions.

% final proporiton completed
conditionMet = (finalprop(:,1) -finalprop(:,2))<=0.001;
% Find the indices of rows where the condition is violated
rowsViolating = find(~conditionMet);
% result: YES!

% final reward
conditionMet = (U_task_with(:,1) - U_task_with(:,2))<0.001;

rowsViolating = find(~conditionMet);
% result: YES!

% total cost
conditionMet = (cost_with(:,1) - cost_with(:,2))<0.001;

rowsViolating = find(~conditionMet);
% result: NO!

% net utility
conditionMet = (Utotal_with(:,1) - Utotal_with(:,2))<0.001;

rowsViolating = find(~conditionMet);
% result: NO!


%% percentage of agents start to work before the first intermediate deadline

sum(Ndaysofprocras(:,1)<7)/length(Ndaysofprocras(:,1)) % 0.1726

sum(Ndaysofprocras(:,2)<7)/length(Ndaysofprocras(:,2)) % 0.9187

nanmean(Ndaysofprocras,1) % 13.6292    4.7328
%% mean of all the outcomes

figure

subplot(1,4,1)
b=bar([1,2],nanmean(Ndaysofprocras,1),'FaceColor','flat');
b.CData(1,:) = [27,158,119]/255;
b.CData(2,:) = [217,95,2]/255;
b.EdgeColor = 'flat';
xticks([1,2])
%xticklabels({'delayed reward','reward upon task completion','reward each milestone (n\in [2,5])','reward each unit of proress'})
xticklabels({'single deadline','interim deadlines'})
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
data = [meanU_task(1),meanCost(1),meanUnet(1);meanU_task(2),meanCost(2),meanUnet(2)];
b = bar(data, 'grouped');

% Set the colors for each bar
for i = 1:size(data, 2)  % Loop over the number of bars in each group
    b(i).FaceColor = 'flat';
    b(i).CData = repmat(threecolor(i, :), size(data, 1), 1);
    b(i).EdgeColor = 'flat';
end
legend({'final reward', 'total cost', 'net utility'})
legend boxoff
xticks([1,2])
xticklabels({'single deadline','interim deadlines'})
yticks(0:0.2:1);
xtickangle(45);
ylabel('utility')
ylim([0,1])
box off
set(gca,'TickDir','out');

subplot(1,4,2)
b=bar([1,2],nanmean(finisheday,1),'FaceColor','flat');
b.CData(1,:) = [27,158,119]/255;
b.CData(2,:) = [217,95,2]/255;
b.EdgeColor = 'flat';
xticks([1,2])
xticklabels({'single deadline','interim deadlines'})
yticks(0:2:TAll(1));
xtickangle(45);
ylabel('average task completion day')
ylim([0,TAll(1)])
box off
set(gca,'TickDir','out');

subplot(1,4,3)
b=bar([1,2],nanmean(finalprop,1),'FaceColor','flat');
b.CData(1,:) = [27,158,119]/255;
b.CData(2,:) = [217,95,2]/255;
b.EdgeColor = 'flat';
xticks([1,2])
xticklabels({'single deadline','interim deadlines'})
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

%% all the outcome histogram
figure
subplot(1,3,1)
histogram(Ndaysofprocras(:,1),'Normalization', 'probability','FaceColor', [27,158,119]/255);
hold on; 
histogram(Ndaysofprocras(:,2),'Normalization', 'probability','FaceColor',[217,95,2]/255);
hold off
legend('single deadline','interim deadlines')
legend boxoff
xticks([0,7,14,21]);
ylim([0,1])
xlabel('days of procrastination')
ylabel('percentage of total (%)')
box off
set(gca,'TickDir','out');

subplot(1,3,2)
histogram(finisheday(:,1),'Normalization', 'probability','FaceColor', [27,158,119]/255);
hold on; 
histogram(finisheday(:,2),'Normalization', 'probability','FaceColor',[217,95,2]/255);
hold off
legend('single deadline','interim deadlines')
legend boxoff
xticks([7,14,21]);
ylim([0,1])
xlabel('task completion day')
ylabel('percentage of total (%)')
box off
set(gca,'TickDir','out');

subplot(1,3,3)
histogram(finalprop(:,1),'Normalization', 'probability','FaceColor', [27,158,119]/255);
hold on; 
histogram(finalprop(:,2),'Normalization', 'probability','FaceColor',[217,95,2]/255);
hold off
legend('single deadline','interim deadlines')
legend boxoff
yticks(0:0.2:1);
ylim([0,1])
xlabel('final proportion completed')
ylabel('percentage of total (%)')
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=900;
height=250;
set(gcf,'position',[x0,y0,width,height])


