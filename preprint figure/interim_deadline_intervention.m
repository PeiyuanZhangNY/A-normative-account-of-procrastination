
close all;
clear 
rng('shuffle')

%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

tic
%% time course of progress
% parameters specific related to two research articles

T = 21; % 21 days as an example for Ariely 2022.
%T = 15; % 15 weeks for Wesp 1986.
N_milestone = 3; % as an example for Ariely 2022
%N_milestone = T;% Wesp1986

deltas=0.01;
k=1;
NMonte = 10000;

alphaAll = rand(1,NMonte); % fixed at 1 for Ariely 2022
%Tmax=14;
%TAll = randi([6,Tmax],1,NMonte));
TAll = T*ones(1,NMonte);  
gammaAll = rand(1,NMonte);
loglambdaAll = rand(1,NMonte); 
lambdaAll = 10.^(loglambdaAll); % lambda>1
%betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
%betaAll = 10.^betaAlllog;
betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
betaAll = 10.^betaAlllog;% fixed at 1 for Ariely 2022
c1All = 2*rand(1,NMonte);
JAll = 0.01*alphaAll(1,NMonte).*rand(1,NMonte);
%N_milestone = randi([2,5],1,NMonte); % has to be smaller than TAll

Npara= 2;% delayed reward, interim deadline 
Ndaysofprocras = nan(NMonte,Npara);
MUCD = nan(NMonte,Npara); % mean unit completion day
%Normalizeddaysofprocras = nan(NMonte,Npara);
%Nnetutility = nan(NMonte, Npara);
finisheday = nan(NMonte,Npara);
finalprop = nan(NMonte,Npara);

WithoptActSeqMatrix1 = nan(NMonte,TAll(1)); % delayed reward
WithoptActSeqMatrix2 = nan(NMonte,TAll(1)); % interim deadline

net_earning_all = nan(NMonte,1);

penalty_unit_all = 0.1*rand(1,NMonte); % 0.03 is for Ariely and Wertenbroch 2002
for nMonte = 1: NMonte
    nMonte
    
    [temp1,] = OptActStateSeq( deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
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
    
    
    
    [temp2,~,net_earning] = OptActStateSeq_interimdeadline( deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),N_milestone,penalty_unit_all(nMonte)); 
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

toc % NMonte=1000; 7 minutes

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

% alternatively, I can calculate net_earning_all directly from the
% following code (for this code, I can also get the total loss due to
% passed interim deadlines)
net_earning_all = nan(NMonte,1);
loss_all = nan(NMonte,1);
reward_interval_list = cumsum(ones(1,N_milestone)*1/N_milestone); %[1/3,2/3,1]
interim_deadline_list = cumsum(ones(1,N_milestone)*round(T/N_milestone)); % [7,14,21]

for nMonte=1:NMonte
    OptActSequence = WithoptActSeqMatrix2(nMonte,:);
    net_earning = 0; % the utility gained in the end - cumulative penalty_loss along the way
    % when t between 1/3*T (7) and 2/3*T (14), if <1/3, then lose penalty_unit
    for t = interim_deadline_list(1):interim_deadline_list(2)-1
        if OptActSequence(t)<1/N_milestone
            net_earning = net_earning -penalty_unit_all(nMonte);
        end
    end

    % when t between 2/3*T (14) and T (21), if <1/3, then lose 2*penalty_unit,
    % if <2/3 but >=1/3, lose penalty_unit
    cumsum_OptActSequence = cumsum(OptActSequence);
    for t = interim_deadline_list(2):T-1
        if cumsum_OptActSequence(t)<1/N_milestone
            net_earning = net_earning -2*penalty_unit_all(nMonte);
        elseif cumsum_OptActSequence(t)<2/N_milestone
            net_earning = net_earning -penalty_unit_all(nMonte);
        end
    end
    % when t=T
    listminusstate = reward_interval_list - cumsum_OptActSequence(T); 
    closestindexOnLeftOfState = sum(listminusstate<=0); %
    net_earning = net_earning-(N_milestone-closestindexOnLeftOfState)*penalty_unit_all(nMonte);
    loss_all(nMonte)=net_earning;
    net_earning_all(nMonte) = net_earning+alphaAll(nMonte)*cumsum_OptActSequence(T)^betaAll(nMonte);
end

% save('interim_deadline_revised.mat') % for single deadline, low level, median
% level and high level penalty_unit
% data is saved save('Ariely2022_1.mat') and save('Ariely2022_2.mat') when N_mileston=3, T=21;
% or save('Wesp1986_1.mat') and save('Wesp1986_2.mat') % when N_milestone=T, T=15;

%save('interim_deadline_revised.mat')
%% time course of progress
load('interim_deadline_revised.mat')

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
% result: YES! ?

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

sum(Ndaysofprocras(:,1)<7)/length(Ndaysofprocras(:,1)) %0.0154

sum(Ndaysofprocras(:,2)<7)/length(Ndaysofprocras(:,2)) %0.7662
%% plots
% using two contradictory color
% example time course of progress
figure
subplot(1,2,1)
%greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
T = TAll(1);
plot(1:T,WithoptActSeqMatrix1(6,:),'o-','Color',[27,158,119]/255);%greencolor(3,:)
hold on
plot(1:T,WithoptActSeqMatrix2(6,:),'o-','Color',[217,95,2]/255);%greencolor(6,:)
hold off
legend('single deadline', 'interim deadline');
legend boxoff

%legend('single deadline (control condition)','interim deadline')
%legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
xticks([7,14,21])
yticks(0:0.2:0.6)
ylim([0,0.6])
box off
set(gca,'TickDir','out');
%axis square

subplot(1,2,2)
T = TAll(1);
plot(1:T,WithoptActSeqMatrix1(183,:),'o-','Color',[27,158,119]/255);
hold on
plot(1:T,WithoptActSeqMatrix2(183,:),'o-','Color',[217,95,2]/255);
hold off


legend('single deadline', 'interim deadline');
legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
xticks([7,14,21])
yticks(0:0.1:0.3)
ylim([0,0.3])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=600;
height=250;
set(gcf,'position',[x0,y0,width,height])

%% patterns of time course of progress
% figure
% for i=1:100
%     if finalprop(i,2)==1
%         plot(1:finisheday(i,2),WithoptActSeqMatrix2(i,1:finisheday(i,2)),'k')
%     else
%         plot(1:21,WithoptActSeqMatrix2(i,:),'k')
%     end
%     hold on
% end
% only two pattern: three ramping activities concatenated and
% each ramping activity is towards each interim deadline.;
% second is one ramping activity toward the final deadline.

% %% when the first pattern of time course of progress, when the second pattern
% % hypothesis: temporal discounting increases, more towards the second
% % pattern
% T = 21; % 21 days for Ariely 2022.
% N_milestone = 3;
% deltas = 0.01;
% k = 1;
% %TAll = randi([Tmin Tmax],1,NMonte);
% alpha = 1;% fixed %rand(1,NMonte);
% gamma = gammaAll(15);%rand(1,NMonte); 
% c1 = c1All(15);%2*rand(1,NMonte);
% %loglambdaAll = rand(1,NMonte); 
% %lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
% lambda = lambdaAll(15); 
% %betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
% %betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
% beta =1;%fixed
% J=JAll(15);
% 
% [temp1,] = OptActStateSeq_Ariely2022_TaskErrorCorrection( deltas,T,k,J,c1,lambda,gamma,alpha,beta,N_milestone);
% [temp2,] = OptActStateSeq_Ariely2022_TaskErrorCorrection( deltas,T,k,J,c1,lambda,1,alpha,beta,N_milestone);
% % temp1: [0,0,0,0,0,0,0.170000000000000,0.170000000000000,0,0,0,0,0,0.160000000000000,0.170000000000000,0,0,0,0,0.160000000000000,0.170000000000000]
% % temp2: [0,0,0,0,0,0,0.200000000000000,0.200000000000000,0,0,0,0,0,0.200000000000000,0.200000000000000,0,0,0,0,0,0.200000000000000]
% % so with increasing gamma, still the same pattern of progress
% [temp2,] = OptActStateSeq_Ariely2022_TaskErrorCorrection( deltas,T,k,J,c1,lambda+6,gamma,alpha,beta,N_milestone);
% % [0,0,0,0,0,0,0.340000000000000,0,0,0,0,0,0,0.330000000000000,0,0,0,0,0,0,0.330000000000000]
% % so with increasing lambda, still the same pattern of progress
% 
% 
% T = 21; % 21 days for Ariely 2022.
% N_milestone = 3;
% deltas = 0.01;
% k = 1;
% %TAll = randi([Tmin Tmax],1,NMonte);
% alpha = 1;% fixed %rand(1,NMonte);
% gamma = gammaAll(3);%rand(1,NMonte);  % 0.127
% c1 = c1All(3);%2*rand(1,NMonte);
% %loglambdaAll = rand(1,NMonte); 
% %lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
% lambda = lambdaAll(3); 
% %betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
% %betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
% beta =1;%fixed
% J=JAll(3);
% 
% [temp1,] = OptActStateSeq_Ariely2022_TaskErrorCorrection( deltas,T,k,J,c1,lambda,gamma,alpha,beta,N_milestone);
% [temp2,] = OptActStateSeq_Ariely2022_TaskErrorCorrection( deltas,T,k,J,c1,lambda,1,alpha,beta,N_milestone);
% % temp1 [0,0,0,0,0,0,0.340000000000000,0,0,0,0,0,0,0.330000000000000,0,0,0,0,0,0,0.330000000000000]
% % temp2 [0,0,0,0,0,0,0.500000000000000,0,0,0,0,0,0,0.500000000000000,0,0,0,0,0,0,0]
% % so with increasing gamma, still the same pattern of progress


%% draw figure
% figure
% subplot(1,4,1)
% scatter(Ndaysofprocras(:,1),Ndaysofprocras(:,2),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:3:TAll(1));
% yticks(0:3:TAll(1));
% xlabel('delayed reward')
% ylabel('interim deadlines')
% axis square
% title('days of procrastination')
% 
% subplot(1,4,3)
% scatter(MUCD(:,1),MUCD(:,2),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:3:TAll(1));
% yticks(0:3:TAll(1));
% xlabel('delayed reward')
% ylabel('interim deadlines')
% axis square
% title('mean unit completion day')
% 
% subplot(1,4,2)
% scatter(finisheday(:,1),finisheday(:,2),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:3:TAll(1));
% yticks(0:3:TAll(1));
% xlabel('delayed reward')
% ylabel('interim deadlines')
% axis square
% title('task completion day')
% 
% subplot(1,4,4)
% scatter(finalprop(:,1),finalprop(:,2),'k')
% hold on
% plot([0,1],[0,1],'k-')
% xlim([0,1])
% ylim([0,1])
% xticks(0:0.2:1);
% yticks(0:0.2:1);
% xlabel('delayed reward')
% %ylabel('Ariely 2022')
% ylabel('interim deadlines')
% axis square
% title('final proportion completed')
% 
% x0=10;
% y0=10;
% width=900;
% height=250;
% set(gcf,'position',[x0,y0,width,height])

%% mean of all the outcomes
%load('interim_deadline_revised.mat')
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
figure
subplot(1,4,1)
b = bar([1,2],nanmean(Ndaysofprocras,1),'FaceColor','flat');
b.CData(1,:) = [27,158,119]/255;%greencolor(3,:);
b.CData(2,:) = [217,95,2]/255;%greencolor(6,:);
b.EdgeColor = 'flat';
xticks([1,2])
%xticklabels({'delayed reward','reward upon task completion','reward each milestone (n\in [2,5])','reward each unit of proress'})
xticklabels({'control','interim deadlines'})
yticks(0:3:TAll(1));
xtickangle(45);
ylabel('average days of procrastination')
ylim([0,TAll(1)])
box off
set(gca,'TickDir','out');

subplot(1,4,4)
% two groups of three bars
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
xticks([1,2,3,4])
xticklabels({'control','interim deadlines'})
yticks(0:0.2:1);
xtickangle(45);
ylabel('utility')
ylim([0,1])
box off
set(gca,'TickDir','out');

subplot(1,4,2)
b=bar([1,2],nanmean(finisheday,1),'FaceColor','flat');
b.CData(1,:) = [27,158,119]/255;%greencolor(3,:);
b.CData(2,:) = [217,95,2]/255;%greencolor(6,:);
b.EdgeColor = 'flat';
xticks([1,2])
xticklabels({'control','interim deadlines'})
yticks(0:3:TAll(1));
xtickangle(45);
ylabel('average task completion day')
ylim([0,TAll(1)])
box off
set(gca,'TickDir','out');

subplot(1,4,3)
b=bar([1,2],nanmean(finalprop,1),'FaceColor','flat');
b.CData(1,:) = [27,158,119]/255;%greencolor(3,:);
b.CData(2,:) = [217,95,2]/255;%greencolor(6,:);
b.EdgeColor = 'flat';
xticks([1,2])
xticklabels({'control','interim deadlines'})
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
histogram(Ndaysofprocras(:,1),'Normalization', 'probability','FaceColor', [27,158,119]/255,'EdgeColor', [27, 158, 119]/255);
hold on; 
histogram(Ndaysofprocras(:,2),'Normalization', 'probability','FaceColor',[217,95,2]/255,'EdgeColor',[217,95,2]/255);
hold off
legend('control','interim deadlines')
legend boxoff
xticks([0,7,14,21]);
ylim([0,1])
xlabel('days of procrastination')
ylabel('percentage of total')
box off
set(gca,'TickDir','out');

subplot(1,3,2)
% Remove NaNs
data1 = finisheday(:,1);
data1 = data1(~isnan(data1));
data2 = finisheday(:,2);
data2 = data2(~isnan(data2));

% Define consistent bin edges for comparison
min_val = min([data1; data2]);
max_val = max([data1; data2]);
binEdges = min_val:1:max_val+1;  % adjust bin width as needed

% Compute histogram counts
[counts1, ~] = histcounts(data1, binEdges);
[counts2, ~] = histcounts(data2, binEdges);

% Normalize to probability
prob1 = counts1 / sum(counts1);
prob2 = counts2 / sum(counts2);

% Plot as bar charts (use bin centers)
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

bar(binCenters, prob1, 'FaceColor', [27,158,119]/255, 'EdgeColor', [27,158,119]/255, 'BarWidth', 1);
hold on;
bar(binCenters, prob2, 'FaceColor', [217,95,2]/255, 'EdgeColor', [217,95,2]/255, 'BarWidth', 0.6);
hold off;

%legend('single deadline','interim deadlines')
%legend boxoff
xticks([7,14,21]);
ylim([0,1])
xlabel('task completion day')
ylabel('percentage of total')
box off
set(gca,'TickDir','out');

subplot(1,3,3)
histogram(finalprop(:,1),'Normalization', 'probability','FaceColor', [27,158,119]/255,'EdgeColor', [27, 158, 119]/255);
hold on; 
histogram(finalprop(:,2),'Normalization', 'probability','FaceColor',[217,95,2]/255,'EdgeColor',[217,95,2]/255);
hold off
%legend('single deadline','interim deadlines')
%legend boxoff
ylim([0,1])
xlabel('final proportion completed')
ylabel('percentage of total')
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=900;
height=250;
set(gcf,'position',[x0,y0,width,height])

