close all;
clear
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

%%%%%%%% first run DelayedRewardAlpha.m code 
%%%%%%% parameter set as follows: deltas = 0.01; k=1;u=0;
%alphaVec=0:0.1:1; lambda=3; beta=2;c1 = 10; gamma=0.4;T=10;

load('alpha.mat')

figure
%% plot two examples of time course curve of progress 
subplot(1,4,1)
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

plot(1:T,WithoptActSeqMatrix(:,3),'o-','Color',greencolor(2,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,6),'o-','Color',greencolor(4,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,9),'o-','Color',greencolor(6,:));
hold off

legend('\alpha=0.2','\alpha=0.5','\alpha=0.8','Interpreter','tex')
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
legend boxoff
xlim([1,T])
xticks(1:T)
box off
set(gca,'TickDir','out');
%% days of procrastination
subplot(1,4,2)
numberzero = nan(1,length(alphaVec));
for alphaIdx = 1:length(alphaVec)
    numberzero(alphaIdx)= nnz(~WithoptActSeqMatrix(:,alphaIdx));
end
plot(alphaVec,numberzero,'ko-')
xlabel('maximum task reward \alpha','Interpreter','tex')
ylabel('days of procrastination')
ylim([0 T])
set(gca, 'XTick', alphaVec(1:2:end));
xlim([alphaVec(1),alphaVec(end)])
box off
set(gca,'TickDir','out');
%% final proportional completed
subplot(1,4,3)
plot(alphaVec,finalprop,'ko-');
xlabel('maximum task reward \alpha','Interpreter','tex')
ylabel('final proportion completed')
ylim([0 1])
set(gca, 'XTick', alphaVec(1:2:end));
xlim([alphaVec(1),alphaVec(end)])
box off
set(gca,'TickDir','out');
%% plot the total cost
subplot(1,4,4)
threecolor = [27,158,119;217,95,2;117,112,179]/255;
plot(alphaVec,U_task_with,'o-','Color',threecolor(1,:))
hold on
plot(alphaVec,cost_with,'o-','Color',threecolor(2,:));
hold on
plot(alphaVec,Utotal_with,'o-','Color',threecolor(3,:))
legend('final reward','total cost','net utility')
legend boxoff
xlabel('maximum task reward \alpha','Interpreter','tex')
ylabel('utility')
ylim([0 1])
set(gca, 'XTick', alphaVec(1:2:end));
xlim([alphaVec(1),alphaVec(end)])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])
% %% plot the U_net
% subplot(1,4,4)
% plot(alphaVec,Utotal_with,'ko-');
% 
% xlabel('maximum task reward \alpha','Interpreter','tex')
% 
% ylabel('net utility')
% ylim([0 1])
% set(gca, 'XTick', alphaVec(1:2:end));
% xlim([alphaVec(1),alphaVec(end)])