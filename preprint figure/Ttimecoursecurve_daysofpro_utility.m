close all;
clear
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

%%%%%%%% first run DelayedRewardT.m code  

load('T.mat')

greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 

% tempplot5 = [0.001,WithoptActSeqMatrix(2:end,8)'];
% tempplot6 = [0.001,0.003, WithoptActSeqMatrix(3:end,9)'];

figure
subplot(1,4,1)
plot(1:10,WithoptActSeqMatrix(1:10,2),'o-','Color',greencolor(2,:));
hold on
%plot(1:T,WithoptActSeqMatrix(:,1),'o--','Color',redcolor(2,:));
%hold on
plot(1:10,WithoptActSeqMatrix(1:10,3),'o-','Color',greencolor(4,:));
hold on
plot(1:10,WithoptActSeqMatrix(1:10,4),'o-','Color',greencolor(6,:));
hold off
legend('T=3','T=4','T=5','Interpreter','tex')
legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')

xlim([1,10])
xticks(1:2:TVec(end))
ylim([0 0.3])
set(gca, 'XTick', 1:10);
box off
set(gca,'TickDir','out');
% figure
% %% plot two examples of time course curve of progress 
% % plot(1:T,tempplot1,'o-','Color',greencolor(1,:));
% % hold on
% % plot(1:T,tempplot2,'o-','Color',greencolor(2,:));
% % hold on
% % plot(1:T,tempplot3,'o-','Color',greencolor(3,:));
% % hold on
% plot(1:T,WithoptActSeqMatrix(:,14),'o-','Color',greencolor(2,:));
% hold on
% plot(1:T,WithoptActSeqMatrix(:,15),'o-','Color',greencolor(4,:));
% hold on
% plot(1:T,WithoptActSeqMatrix(:,16),'o-','Color',greencolor(6,:));
% hold off
% 
% legend('T=15','T=16','T=17','Interpreter','tex')
% xlabel('time \itt','Interpreter','tex')
% ylabel('effort \ita','Interpreter','tex')
% xlim([1,TVec(end)])
% % xticks(1:2:TVec(end))
% xticks(1:TVec(end)) % for zoom
% % ylim([0 0.3])
% ylim([0 0.05]) % for zoom
% yticks(0:0.01:0.05) % for zoom

%% normalized days of procrastination
subplot(1,4,2)
numberzero = nan(1,length(TVec));
for TIdx = 1:9
    numberzero(TIdx)= nnz(~WithoptActSeqMatrix(1:TVec(TIdx),TIdx))/TVec(TIdx);
end
plot(TVec,numberzero,'ko-')
xlabel('total time \itT','Interpreter','tex')
ylabel('days of procrastination/\itT','Interpreter','tex')
ylim([0 1])
xlim([2,10])
xticks(2:1:10)
box off
set(gca,'TickDir','out');

%% final proportional completed
subplot(1,4,3)
plot(TVec,finalprop,'ko-');
xlabel('total time \itT','Interpreter','tex')
ylabel('final proportion completed')
ylim([0 1])
xlim([2,10])
xticks(2:1:10)
box off
set(gca,'TickDir','out');

%% plot the total cost
subplot(1,4,4)
threecolor = [27,158,119;217,95,2;117,112,179]/255;
plot(TVec,U_task_with,'o-','Color',threecolor(1,:))
hold on
plot(TVec,cost_with,'o-','Color',threecolor(2,:))
hold on
plot(TVec,Utotal_with,'o-','Color',threecolor(3,:));
legend('final reward','total cost','net utility')
legend boxoff
xlabel('total time \itT')
ylabel('utility')
ylim([0 1])
xlim([2,10])
xticks(2:1:10)
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])

