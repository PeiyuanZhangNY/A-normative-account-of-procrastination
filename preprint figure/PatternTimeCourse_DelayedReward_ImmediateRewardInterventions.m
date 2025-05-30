%%%%%%%%%%% the agent works little in the early days, and work more and more towards the deadline.
close all;
clear

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

%% simulate the time course
%% the following parameter is used to plot PatternTimeCourse_DelayedReward_ImmediateRewardInterventions.pdf
deltas = 0.01;
k = 1;
T=10;
%TAll = randi([Tmin Tmax],1,NMonte);
alpha = 0.5;%rand(1,NMonte);
gamma = 0.3;%rand(1,NMonte); 
c1 = 2;%2*rand(1,NMonte);
%loglambdaAll = rand(1,NMonte); 
%lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
lambda = 3; 
%betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
%betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
beta =1;
J=0;
%JAll = zeros(1,NMonte);
%JAll = 0.01*alphaAll.*rand(1,NMonte);

%% the following parameter is used to plot PatternTimeCourse_DelayedReward_ImmediateRewardInterventions_2.pdf
% deltas = 0.01;
% k = 1;
% T=10;
% %TAll = randi([Tmin Tmax],1,NMonte);
% alpha = 0.5;%rand(1,NMonte);
% gamma = 0.4;%rand(1,NMonte); 
% c1 = 2;%2*rand(1,NMonte);
% %loglambdaAll = rand(1,NMonte); 
% %lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
% lambda = 3; 
% %betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
% %betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
% beta =2;
% J=0;
% %JAll = zeros(1,NMonte);
% %JAll = 0.01*alphaAll.*rand(1,NMonte);

%% the following parameter is used to plot PatternTimeCourse_DelayedReward_ImmediateRewardInterventions_3.pdf
% deltas = 0.01;
% k = 1;
% T=10;
% %TAll = randi([Tmin Tmax],1,NMonte);
% alpha = 0.8;%rand(1,NMonte);
% gamma = 0.4;%rand(1,NMonte); 
% c1 = 2;%2*rand(1,NMonte);
% %loglambdaAll = rand(1,NMonte); 
% %lambdaAll = 10.^(loglambdaAll);% lambda>1, 1<lambda<10
% lambda = 3; 
% %betaAlllog = log10(0.1)+(log10(10)-log10(0.1))*rand(1,NMonte);
% %betaAll = 10.^betaAlllog; % 0.1<beta<10 so beta can be <lambda, or >lambda
% beta =2;
% J=0;

% delayed reward

[DelayedRewardTimeCourse, ]=OptActStateSeq(deltas,T,k,J,c1,lambda,gamma,alpha,beta); 


% RewardCompletionTime

%for nMonte = 1:NMonte
    %nMonte
[RewardCompletionTimeTimeCourse, ]=OptActStateSeqRewardCompletionTime(deltas,T,k,J,c1,lambda,gamma,alpha,beta); 
if ~isempty(find(RewardCompletionTimeTimeCourse==0,1))
    RewardCompletionTimeCompletionDay = find(RewardCompletionTimeTimeCourse==0,1)-1;
else
    RewardCompletionTimeCompletionDay=T;
end
%end

% Reward each Milestone
N_milestone = 4;
%for nMonte = 1:NMonte
    %nMonte
[RewardMilestoneTimeCourse, ]=OptActStateSeqRewardMilestone(deltas,T,k,J,c1,lambda,gamma,alpha,beta,N_milestone); 
if ~isempty(find(RewardMilestoneTimeCourse==0,1))
    RewardMilestoneCompletionDay = find(RewardMilestoneTimeCourse==0,1)-1;
else
    RewardMilestoneCompletionDay=T;
end
%end

% Reward each unit of progress

%for nMonte = 1:NMonte
    %nMonte
[RewardUnitofProgressTimeCourse, ]=OptActStateSeqRewardEachUnitOfProgress(deltas,T,k,J,c1,lambda,gamma,alpha,beta); 
if ~isempty(find(RewardUnitofProgressTimeCourse==0,1))
    RewardUnitofProgressCompletionDay = find(RewardUnitofProgressTimeCourse==0,1)-1;
else
    RewardUnitofProgressCompletionDay=T;
end
%end

CombinedTimeCourseFromAllconditions = [DelayedRewardTimeCourse,RewardCompletionTimeTimeCourse,RewardMilestoneTimeCourse,RewardUnitofProgressTimeCourse];


%% draw figure to show the time course  

% time course of progress
figure
%subplot(1,3,1)
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

plot(1:RewardUnitofProgressCompletionDay,RewardUnitofProgressTimeCourse(1:RewardUnitofProgressCompletionDay),'o-','Color',greencolor(2,:));
hold on
plot(1:RewardMilestoneCompletionDay,RewardMilestoneTimeCourse(1:RewardMilestoneCompletionDay),'o-','Color',greencolor(4,:));
hold on
plot(1:RewardCompletionTimeCompletionDay,RewardCompletionTimeTimeCourse(1:RewardCompletionTimeCompletionDay),'o-','Color',greencolor(6,:));
hold on
plot(1:T,DelayedRewardTimeCourse,'ko-');
hold off

legend('reward each unit of progress','reward each milestone','reward upon task completion','delayed reward')
legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
xticks(1:T)
yticks(0:0.1:0.5)
box off
set(gca,'TickDir','out');
%axis square

x0=10;
y0=10;
width=300;
height=250;
set(gcf,'position',[x0,y0,width,height])

% time course of cumulative progress
% figure
% greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% 
% plot(1:T,cumsum(RewardUnitofProgressTimeCourse),'o-','Color',greencolor(2,:));
% hold on
% plot(1:T,cumsum(RewardMilestoneTimeCourse),'o-','Color',greencolor(4,:));
% hold on
% plot(1:T,cumsum(RewardCompletionTimeTimeCourse),'o-','Color',greencolor(6,:));
% hold on
% plot(1:T,cumsum(DelayedRewardTimeCourse),'ko-');
% hold off
% 
% %legend('reward each unit of progress','reward each milestone','reward upon task completion','delayed reward')
% legend boxoff
% xlabel('time \itt','Interpreter','tex')
% ylabel('cumulative progress \its','Interpreter','tex')
% xlim([1,T])
% xticks(1:T)
% yticks(0:0.2:1)


% %% test whether the conclusion (work more and more towards the deadline) hold true
% tplus1minust = nan(1,NMonte*Tmax);
% n=1; 
% for nMonte = 1:NMonte
%      effortvec = WithoptActSeqMatrix(nMonte,1:TAll(nMonte));
%      for effortidx = 1:length(effortvec)-1
%          tplus1minust(1,n)=effortvec(effortidx+1)-effortvec(effortidx);
%          if tplus1minust(1,n)<0
%              disp(nMonte)
%          end
%              
%          n=n+1;
%      end
% end
% tplus1minust(tplus1minust<0)
% figure
% histogram(tplus1minust)

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

