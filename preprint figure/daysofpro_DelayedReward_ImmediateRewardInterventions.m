%%%%%%%%% 
% this is revised from
% daysofpro_utility_DelayedReward_RewardCompetionTime.m
close all;
clear 
rng('shuffle')
% data saved in 'ImmediateRewardIntervention_2.mat'
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);
tic
%% time course of progress

deltas=0.01;
k=1;
NMonte = 10000;

alphaAll = rand(1,NMonte);
%Tmax=14;
%TAll = randi([6,Tmax],1,NMonte);
TAll = 10*ones(1,NMonte); %fixed at 10 days. 
gammaAll = rand(1,NMonte);
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
    
    
    
    [temp2,] = OptActStateSeqRewardCompletionTime( deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
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
    
    [temp3,] = OptActStateSeqRewardMilestone( deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte),N_milestone(nMonte)); 
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

    [temp4,] = OptActStateSeqRewardEachUnitOfProgress( deltas,TAll(nMonte),k,JAll(nMonte),c1All(nMonte),lambdaAll(nMonte),gammaAll(nMonte),alphaAll(nMonte),betaAll(nMonte)); 
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

checknMonte = nan(1,100);
i=0;
for nMonte=1:NMonte
    if finalprop(nMonte,1)>finalprop(nMonte,3)
        i=i+1;
        checknMonte(i)=nMonte;
    end
end

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


%% save data
% the result is saved as ImmediateRewardIntervention_2.mat
load('ImmediateRewardIntervention_2.mat')

%% check if the average result hold true for each simulation
% days of prorastination delayed>upon task completion>upon milestones>upon
% unit of work
conditionMet = Ndaysofprocras(:,1) >= Ndaysofprocras(:,2) & ...
               Ndaysofprocras(:,2) >= Ndaysofprocras(:,3) & ...
               Ndaysofprocras(:,3) >= Ndaysofprocras(:,4);

% Find the indices of rows where the condition is violated
rowsViolating = find(~conditionMet);
% result: 0×1 empty; YES!  In terms of days of procrastination, in every single simulation, offering immediate rewards reduces the number of days of delay. The higher the immediacy level of interventions, the fewer days of procrastination.

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
% result: NO!

%% figure: patterns of time course of progress

% % % reward upon task completion
% % figure
% % for i=1:100
% %     if finalprop(i,2)==1
% %         plot(1:finisheday(i,2),WithoptActSeqMatrix2(i,1:finisheday(i,2)),'k')
% %     else
% %         plot(1:10,WithoptActSeqMatrix2(i,:),'k')
% %     end
% %     hold on
% % end
% % % shape: ramping up activity until the task completion day
% % 
% % % reward each milestone
% % figure
% % for i=1:100
% %     if finalprop(i,3)==1
% %         plot(1:finisheday(i,3),WithoptActSeqMatrix3(i,1:finisheday(i,3)),'k')
% %     else
% %         plot(1:10,WithoptActSeqMatrix3(i,:),'k')
% %     end
% %     hold on
% %     pause
% % end
% % 
% % % example time course of progress

% % figure for (patterns of time course of progress in reward-each-milestone schedule 
% set(0,'DefaultLineLineWidth',2);
% set(0,'DefaultAxesFontSize',16);
% figure
% subplot(2,3,1)
% T=10;
% working steadily
% i=8;
% plotVector = NaN(1, T); 
% plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
% plot(1:T,plotVector,'ko-');
% xlabel('time \itt','Interpreter','tex')
% ylabel('progress \Delta\its','Interpreter','tex')
% xlim([1,T])
% ylim([0,0.8])
% xticks(1:T)
% yticks(0:0.2:0.8)
% set(gca,'TickDir','out');
% axis square
% 
% subplot(2,3,2)
% ramping up 
% i=14;
% plotVector = NaN(1, T); 
% plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
% plot(1:T,plotVector,'ko-');  
% xlabel('time \itt','Interpreter','tex')
% ylabel('progress \Delta\its','Interpreter','tex')
% xlim([1,T])
% ylim([0,0.8])
% xticks(1:T)
% yticks(0:0.2:0.8)
% set(gca,'TickDir','out');
% axis square
% 
% subplot(2,3,3)
% going down
% i = 7;
% plotVector = NaN(1, T); 
% plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
% plot(1:T,plotVector,'ko-');  
% xlabel('time \itt','Interpreter','tex')
% ylabel('progress \Delta\its','Interpreter','tex')
% xlim([1,T])
% ylim([0,0.8])
% xticks(1:T)
% yticks(0:0.2:0.8)
% set(gca,'TickDir','out');
% axis square
% 
% subplot(2,3,4)
% first going down and then ramping up (U=shape)
% i = 49;
% plotVector = NaN(1, T); 
% plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
% plot(1:T,plotVector,'ko-');  
% xlabel('time \itt','Interpreter','tex')
% ylabel('progress \Delta\its','Interpreter','tex')
% xlim([1,T])
% ylim([0,0.8])
% xticks(1:T)
% yticks(0:0.2:0.8)
% set(gca,'TickDir','out');
% axis square
% 
% subplot(2,3,5)
% first ramping up and then going down (inverted U=shape)
% i = 74;
% plotVector = NaN(1, T); 
% plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
% plot(1:T,plotVector,'ko-'); 
% xlabel('time \itt','Interpreter','tex')
% ylabel('progress \Delta\its','Interpreter','tex')
% xlim([1,T])
% ylim([0,0.8])
% xticks(1:T)
% yticks(0:0.2:0.8)
% set(gca,'TickDir','out');
% axis square
% 
% subplot(2,3,6)
% fluctuate
% i = 67;
% plotVector = NaN(1, T); 
% plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
% plot(1:T,plotVector,'ko-'); 
% xlabel('time \itt','Interpreter','tex')
% ylabel('progress \Delta\its','Interpreter','tex')
% xlim([1,T])
% ylim([0,0.8])
% xticks(1:T)
% yticks(0:0.2:0.8)
% set(gca,'TickDir','out');
% axis square


% % 
% % 
% % % reward each unit of progress
% % figure
% % for i=1:100
% %     if finalprop(i,4)==1
% %         plot(1:finisheday(i,4),WithoptActSeqMatrix4(i,1:finisheday(i,4)),'k')
% %     else
% %         plot(1:10,WithoptActSeqMatrix4(i,:),'k')
% %     end
% %     hold on
% %     pause
% % end
% % % decreasing, ramping up, inverted U-shape (by observing first 100
% samples)

% % double check types of patterns of time course of progress with coding
% % check if there is U-shape (strict version and relax version including
% % constant change); Result: we did NOT find U-shape.
% % strict vesion

% for i=1:NMonte
%     if isnan(finisheday(i,4)) 
%         vec = WithoptActSeqMatrix4(i,:);
%     else        
%         vec = WithoptActSeqMatrix4(i,1:finisheday(i,4));
%     end
%     
%     if length(vec)>=3
%         % Find the minimum value and its index
%         [min_val, min_idx] = min(vec);
%         
%         % min cannot be the last or the first, then it would only increase or
%         % decrease
%         
%         if min_idx~=1 && min_idx~=length(vec)
% 
%         % Check for decreasing trend before the minimum
%         before_min = vec(1:min_idx);
%         is_decreasing = all(diff(before_min) < 0);
% 
%         % Check for increasing trend after the minimum
%         after_min = vec(min_idx:end);
%         is_increasing = all(diff(after_min) > 0);
% 
%         % Combine both checks
%         is_U_shape = is_decreasing && is_increasing;
%         if is_U_shape
%             disp(i);
%             disp(vec);
%         end
%         
%         end
%     end
%         
% end

% % relax version
% for i=1:NMonte
%     if isnan(finisheday(i,4)) 
%         vec = WithoptActSeqMatrix4(i,:);
%     else        
%         vec = WithoptActSeqMatrix4(i,1:finisheday(i,4));
%     end
%     
%     if length(vec)>=3
%         % Find the minimum value and its index
%         [min_val, min_idx] = min(vec);
%         
%         % min cannot be the last or the first, then it would only increase or
%         % decrease
%         
%         if min_idx~=1 && min_idx~=length(vec)
% 
%         % Check for decreasing trend before the minimum
%         before_min = vec(1:min_idx);
%         is_decreasing = all(diff(before_min) <= 0) && sum(diff(before_min)<0)>0;
% 
%         % Check for increasing trend after the minimum
%         after_min = vec(min_idx:end);
%         is_increasing = all(diff(after_min) >= 0) && sum(diff(after_min)>0)>0;
% 
%         % Combine both checks
%         is_U_shape = is_decreasing && is_increasing;
%         if is_U_shape
%             disp(i);
%             disp(vec);
%         end
%         
%         end
%     end
%         
% end

% double check types of patterns of time course of progress with coding
% check if there is fluctuate; Result: no fluctuation.

% NaNWithoptActSeqMatrix4 = nan(NMonte,10);
% diffWithoptActSeqMatrix4 = nan(NMonte,9);
% for i=1:NMonte
%     if isnan(finisheday(i,4)) 
%         NaNWithoptActSeqMatrix4(i,:) = WithoptActSeqMatrix4(i,:);
%         diffWithoptActSeqMatrix4(i,:) = diff(WithoptActSeqMatrix4(i,:));
%     else        
%         NaNWithoptActSeqMatrix4(i,1:finisheday(i,4)) = WithoptActSeqMatrix4(i,1:finisheday(i,4));
%         diffWithoptActSeqMatrix4(i,:) = diff(NaNWithoptActSeqMatrix4(i,:));
%     end     
% end


 
% 
% %% figure: compare across all the conditions
% % days of procrastination:
% figure
% subplot(2,3,1) % first row: comparing each immediate reward with delayed reward; 
% % second row: compare acorss immediate reward interventions
% scatter(Ndaysofprocras(:,1),Ndaysofprocras(:,2),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward upon task completion')
% axis square
% box off
% set(gca,'TickDir','out');
% 
% subplot(2,3,2) 
% scatter(Ndaysofprocras(:,1),Ndaysofprocras(:,3),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward each milestone')
% title('days of procrastination')
% axis square
% box off
% set(gca,'TickDir','out');
% 
% subplot(2,3,3) 
% scatter(Ndaysofprocras(:,1),Ndaysofprocras(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward each unit of progress')
% axis square
% box off
% set(gca,'TickDir','out');
% 
% subplot(2,3,4) 
% scatter(Ndaysofprocras(:,2),Ndaysofprocras(:,3),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward upon task completion')
% ylabel('reward each milestone')
% axis square
% box off
% set(gca,'TickDir','out');
% 
% subplot(2,3,5) 
% scatter(Ndaysofprocras(:,2),Ndaysofprocras(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward upon task completion')
% ylabel('reward each unit of progress')
% axis square
% 
% subplot(2,3,6) 
% scatter(Ndaysofprocras(:,3),Ndaysofprocras(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward each milestone')
% ylabel('reward each unit of progress')
% axis square
% 
% x0=10;
% y0=10;
% width=900;
% height=500;
% set(gcf,'position',[x0,y0,width,height])
% 
% % MUCD:
% figure
% subplot(2,3,1) % first row: comparing each immediate reward with delayed reward; 
% % second row: compare acorss immediate reward interventions
% scatter(MUCD(:,1),MUCD(:,2),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward upon task completion')
% axis square
% 
% subplot(2,3,2) 
% scatter(MUCD(:,1),MUCD(:,3),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward each milestone')
% title('mean unit completion day')
% axis square
% 
% subplot(2,3,3) 
% scatter(MUCD(:,1),MUCD(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward each unit of progress')
% axis square
% 
% subplot(2,3,4) 
% scatter(MUCD(:,2),MUCD(:,3),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward upon task completion')
% ylabel('reward each milestone')
% axis square
% 
% subplot(2,3,5) 
% scatter(MUCD(:,2),MUCD(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward upon task completion')
% ylabel('reward each unit of progress')
% axis square
% 
% subplot(2,3,6) 
% scatter(MUCD(:,3),MUCD(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward each milestone')
% ylabel('reward each unit of progress')
% axis square
% 
% x0=10;
% y0=10;
% width=900;
% height=500;
% set(gcf,'position',[x0,y0,width,height])
% 
% % task completion day:
% figure
% subplot(2,3,1) % first row: comparing each immediate reward with delayed reward; 
% % second row: compare acorss immediate reward interventions
% scatter(finisheday(:,1),finisheday(:,2),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward upon task completion')
% axis square
% 
% subplot(2,3,2) 
% scatter(finisheday(:,1),finisheday(:,3),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward each milestone')
% title('task completion day')
% axis square
% 
% subplot(2,3,3) 
% scatter(finisheday(:,1),finisheday(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('delayed reward')
% ylabel('reward each unit of progress')
% axis square
% 
% subplot(2,3,4) 
% scatter(finisheday(:,2),finisheday(:,3),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward upon task completion')
% ylabel('reward each milestone')
% axis square
% 
% subplot(2,3,5) 
% scatter(finisheday(:,2),finisheday(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward upon task completion')
% ylabel('reward each unit of progress')
% axis square
% 
% subplot(2,3,6) 
% scatter(finisheday(:,3),finisheday(:,4),'k')
% hold on
% plot([0,TAll(1)],[0,TAll(1)],'k-')
% xlim([0,TAll(1)])
% ylim([0,TAll(1)])
% xticks(0:2:TAll(1));
% yticks(0:2:TAll(1));
% xlabel('reward each milestone')
% ylabel('reward each unit of progress')
% axis square
% 
% x0=10;
% y0=10;
% width=900;
% height=500;
% set(gcf,'position',[x0,y0,width,height])
% 
% % final proportion completed:
% figure
% subplot(2,3,1) % first row: comparing each immediate reward with delayed reward; 
% % second row: compare acorss immediate reward interventions
% scatter(finalprop(:,1),finalprop(:,2),'k')
% hold on
% plot([0,1],[0,1],'k-')
% xlim([0,1])
% ylim([0,1])
% xticks(0:0.2:1);
% yticks(0:0.2:1);
% xlabel('delayed reward')
% ylabel('reward upon task completion')
% axis square
% 
% subplot(2,3,2) 
% scatter(finalprop(:,1),finalprop(:,3),'k')
% hold on
% plot([0,1],[0,1],'k-')
% xlim([0,1])
% ylim([0,1])
% xticks(0:0.2:1);
% yticks(0:0.2:1);
% xlabel('delayed reward')
% ylabel('reward each milestone')
% title('final proportion completed')
% axis square
% 
% subplot(2,3,3) 
% scatter(finalprop(:,1),finalprop(:,4),'k')
% hold on
% plot([0,1],[0,1],'k-')
% xlim([0,1])
% ylim([0,1])
% xticks(0:0.2:1);
% yticks(0:0.2:1);
% xlabel('delayed reward')
% ylabel('reward each unit of progress')
% axis square
% 
% subplot(2,3,4) 
% scatter(finalprop(:,2),finalprop(:,3),'k')
% hold on
% plot([0,1],[0,1],'k-')
% xlim([0,1])
% ylim([0,1])
% xticks(0:0.2:1);
% yticks(0:0.2:1);
% xlabel('reward upon task completion')
% ylabel('reward each milestone')
% axis square
% 
% subplot(2,3,5) 
% scatter(finalprop(:,2),finalprop(:,4),'k')
% hold on
% plot([0,1],[0,1],'k-')
% xlim([0,1])
% ylim([0,1])
% xticks(0:0.2:1);
% yticks(0:0.2:1);
% xlabel('reward upon task completion')
% ylabel('reward each unit of progress')
% axis square
% 
% subplot(2,3,6) 
% scatter(finalprop(:,3),finalprop(:,4),'k')
% hold on
% plot([0,1],[0,1],'k-')
% xlim([0,1])
% ylim([0,1])
% xticks(0:0.2:1);
% yticks(0:0.2:1);
% xlabel('reward each milestone')
% ylabel('reward each unit of progress')
% axis square
% 
% x0=10;
% y0=10;
% width=900;
% height=500;
% set(gcf,'position',[x0,y0,width,height])

%% mean of all the outcomes
figure
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

subplot(1,4,1)
b=bar([1,2,3,4],mean(Ndaysofprocras,1),'FaceColor','flat');
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
ylabel('days of procrastination')
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
ylabel('task completion day')
ylim([0,TAll(1)])
box off
set(gca,'TickDir','out');

subplot(1,4,3)
b=bar([1,2,3,4],mean(finalprop,1),'FaceColor','flat');
b.CData(1,:) = [0,0,0];
b.CData(2,:) = greencolor(6,:);
b.CData(3,:) = greencolor(4,:);
b.CData(4,:) = greencolor(2,:);
b.EdgeColor = 'flat';
xticks([1,2,3,4])
xticklabels({'delayed reward','reward upon task completion','reward each milestone','reward each unit of proress'})
yticks(0:0.2:1);
xtickangle(45);
ylabel('final proportion completed')
ylim([0,1])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=300;
set(gcf,'position',[x0,y0,width,height])

%% alternative figure 
% it is not comparing between each immediate reward interventions and
% delayed reward. The results are:
% days of procrastination: delayed reward>reward upon task
% completion>reward each mileston>reward each unit of progress
% mean unit completion day:
% task completion day: 
% final proportion completed: 

% figure
% subplot(1,4,1)
% for nMonte = 1:NMonte
%     plot([1,3,5,7],Ndaysofprocras(nMonte,:),'ko-')
%     hold on
% end
% ylabel('days of procrastination','Interpreter','tex')
% xlim([0,8])
% xticks([1,3,5,7])
% xticklabels({'delayed reward','reward upon task completion','reward each milestone (n\in [2,5]','reward each unit of proress'})
% yticks(0:2:TAll(1));
% 
% subplot(1,4,2)
% for nMonte = 1:NMonte
%     plot([1,3,5,7],MUCD(nMonte,:),'ko-')
%     hold on
% end
% ylabel('mean unit completion day','Interpreter','tex')
% xlim([0,8])
% xticks([1,3,5,7])
% xticklabels({'delayed reward','reward upon task completion','reward each milestone (n\in [2,5]','reward each unit of proress'})
% yticks(0:2:TAll(1));
% 
% subplot(1,4,3)
% for nMonte = 1:NMonte
%     plot([1,3,5,7],finisheday(nMonte,:),'ko-')
%     hold on
% end
% ylabel('task completion day','Interpreter','tex')
% xlim([0,8])
% xticks([1,3,5,7])
% xticklabels({'delayed reward','reward upon task completion','reward each milestone (n\in [2,5]','reward each unit of proress'})
% yticks(0:2:TAll(1));
% 
% subplot(1,4,4)
% for nMonte = 1:NMonte
%     plot([1,3,5,7],finalprop(nMonte,:),'ko-')
%     hold on
% end
% ylabel('fina proportion completed','Interpreter','tex')
% xlim([0,8])
% ylim([0,1])
% xticks([1,3,5,7])
% xticklabels({'delayed reward','reward upon task completion','reward each milestone (n\in [2,5]','reward each unit of proress'})
% yticks(0:0.2:1);
% 
% 
% 
