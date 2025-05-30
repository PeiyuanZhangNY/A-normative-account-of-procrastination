close all;
clear
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

%%%%%%%% first run DelayedRewardGamma.m and save the result as gamma.mat
% upper panel
%%%%%%% parameter set as follows: deltas = 0.01; k=1;u=0;eta=7; lambda=3;alpha=1; beta=1;deltagamma=0.2;gammaVec = 0:deltagamma:1;T = 10;
load('gamma.mat')
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

figure
%% plot three examples of time course curve of progress in different gamma
%subplot(1,3,1)
subplot(1,4,1)
plot(1:T,WithoptActSeqMatrix(:,1),'o-','Color',greencolor(2,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,2),'o-','Color',greencolor(4,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,5),'o-','Color',greencolor(6,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,end),'ko-');
hold off

legend('\gamma=0','\gamma=0.1','\gamma=0.4','\gamma=1')
legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
xticks(1:T)
box off
set(gca,'TickDir','out');

%% days of procrastination over gamma  
%subplot(1,3,2)
subplot(1,4,2)
numberzero = nan(1,length(gammaVec));
for gammaIdx = 1:length(gammaVec)
    numberzero(gammaIdx)= nnz(~WithoptActSeqMatrix(:,gammaIdx));
end
plot(gammaVec,numberzero,'ko-')
xlabel('discount rate \it\gamma','Interpreter','tex')
ylabel('days of procrastination')
ylim([0 T])
set(gca, 'XTick', 0:0.2:1);
%set(gca, 'yTick', 0:1:10);
box off
set(gca,'TickDir','out');

%% final proportion of completion  
%subplot(1,3,3)
subplot(1,4,3)
plot(gammaVec,U_task_with,'ko-');
xlabel('discount rate \it\gamma','Interpreter','tex')
ylabel('final proportion completed')
ylim([0 1])
set(gca, 'XTick', 0:0.2:1);

x0=10;
y0=10;
width=900;
height=250;
set(gcf,'position',[x0,y0,width,height])
box off
set(gca,'TickDir','out');

%% total effort cost
subplot(1,4,4)
threecolor = [27,158,119;217,95,2;117,112,179]/255;
plot(gammaVec,U_task_with,'o-','Color',threecolor(1,:))
hold on
plot(gammaVec,cost_with,'o-','Color',threecolor(2,:));
hold on
plot(gammaVec,Utotal_with,'o-','Color',threecolor(3,:))
legend('final reward','total cost','net utility')
legend boxoff
xlabel('discount rate \it\gamma','Interpreter','tex')
ylabel('utility')
ylim([0 1])
set(gca, 'XTick', 0:0.2:1);
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])


% %%%%%%%% second run DelayedRewardGamma.m and save the result as gamma_cost.mat
% % lower pannel
% load('gamma_cost.mat') % eta=0.5, lambda=2; beta=1;
% greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% 
% figure
% %% plot three examples of time course curve of progress in different gamma
% %subplot(1,3,1)
% subplot(1,4,1)
% plot(1:T,WithoptActSeqMatrix(:,1),'o-','Color',greencolor(2,:));
% hold on
% plot(1:T,WithoptActSeqMatrix(:,2),'o-','Color',greencolor(4,:));
% hold on
% plot(1:T,WithoptActSeqMatrix(:,5),'o-','Color',greencolor(6,:));
% hold on
% plot(1:T,WithoptActSeqMatrix(:,end),'ko-');
% hold off
% 
% legend('\gamma=0','\gamma=0.1','\gamma=0.4','\gamma=1')
% legend boxoff
% xlabel('time \itt','Interpreter','tex')
% ylabel('progress \Delta\its','Interpreter','tex')
% xlim([1,T])
% xticks(1:T)
% box off
% %% days of procrastination over gamma  
% %subplot(1,3,2)
% subplot(1,4,2)
% numberzero = nan(1,length(gammaVec));
% for gammaIdx = 1:length(gammaVec)
%     numberzero(gammaIdx)= nnz(~WithoptActSeqMatrix(:,gammaIdx));
% end
% plot(gammaVec,numberzero,'ko-')
% xlabel('discount rate \it\gamma','Interpreter','tex')
% ylabel('days of procrastination')
% ylim([0 T])
% set(gca, 'XTick', 0:0.2:1);
% %set(gca, 'yTick', 0:1:10);
% box off
% 
% %% final proportion of completion  
% %subplot(1,3,3)
% subplot(1,4,3)
% plot(gammaVec,U_task_with,'ko-');
% xlabel('discount rate \it\gamma','Interpreter','tex')
% ylabel('final proportion completed')
% ylim([0 1])
% set(gca, 'XTick', 0:0.2:1);
% 
% x0=10;
% y0=10;
% width=900;
% height=250;
% set(gcf,'position',[x0,y0,width,height])
% box off
% %% total effort cost
% subplot(1,4,4)
% threecolor = [27,158,119;217,95,2;117,112,179]/255;
% plot(gammaVec,U_task_with,'o-','Color',threecolor(1,:))
% hold on
% plot(gammaVec,cost_with,'o-','Color',threecolor(2,:));
% hold on
% plot(gammaVec,Utotal_with,'o-','Color',threecolor(3,:))
% legend('final reward','total cost','net utility')
% legend boxoff
% xlabel('discount rate \it\gamma','Interpreter','tex')
% ylabel('utility')
% ylim([0 0.5])
% set(gca, 'XTick', 0:0.2:1);
% box off
% 
% x0=10;
% y0=10;
% width=1200;
% height=250;
% set(gcf,'position',[x0,y0,width,height])
