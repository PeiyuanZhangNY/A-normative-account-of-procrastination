close all;
clear
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

%%%%%%%% first run DelayedRewardDistractor.m code  
%%%%%%% parameter set as follows: deltas = 0.01; k=1;
%uvec=0:0.002:0.024;eta=7; lambda=3;alpha=1; beta=2;gamma = 0.6; T = 10;

load('J.mat')

figure
%% plot two examples of time course curve of progress 
subplot(1,4,1)
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

plot(1:T,WithoptActSeqMatrix(:,2),'o-','Color',greencolor(2,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,5),'o-','Color',greencolor(4,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,9),'o-','Color',greencolor(6,:));
hold off

legend('J=0.002','J=0.008','J=0.016','Interpreter','tex')
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
legend boxoff
xlim([1,10])
xticks(1:2:T)
ylim([0 0.3])
set(gca, 'XTick', 1:10);
box off
set(gca,'TickDir','out');

%% days of procrastination
subplot(1,4,2)
numberzero = nan(1,length(uvec));
for alphaIdx = 1:length(uvec)
    numberzero(alphaIdx)= nnz(~WithoptActSeqMatrix(:,alphaIdx));
end
plot(uvec,numberzero,'ko-')
xlabel('utility of alternative activities \itJ','Interpreter','tex')
ylabel('days of procrastination')
ylim([0 T])
xticks(uvec(1:4:end))
xlim([uvec(1),uvec(end)])
box off
set(gca,'TickDir','out');

%% final proportional completed
subplot(1,4,3)
plot(uvec,finalprop,'ko-');
xlabel('utility of alternative activities \itJ','Interpreter','tex')
ylabel('final proportion completed')
ylim([0 1])
xlim([uvec(1),uvec(end)])
xticks(uvec(1:4:end))
box off
set(gca,'TickDir','out');

%% plot the total cost
subplot(1,4,4)
threecolor = [27,158,119;217,95,2;117,112,179]/255;
plot(uvec,U_task_with,'o-','Color',threecolor(1,:))
hold on
plot(uvec,cost_with,'o-','Color',threecolor(2,:));
hold on
plot(uvec,Utotal_with,'o-','Color',threecolor(3,:))
legend('final reward','total cost','net utility')
legend boxoff
xlabel('utility of alternative activities \itJ','Interpreter','tex')
ylabel('utility')
ylim([0 1])
set(gca, 'XTick', uvec(1:4:end));
xlim([uvec(1),uvec(end)])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])

% %% plot the U_net
% subplot(1,4,4)
% plot(uvec,Utotal_with,'ko-');
% 
% xlabel('utility of alternative activities','Interpreter','tex')
% 
% ylabel('net utility')
% ylim([0 1])
% set(gca, 'XTick', uvec(1:4:end));
% xlim([uvec(1),uvec(end)])
