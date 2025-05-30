close all;
clear
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

%%%%%%%% first run DelayedRewardCmax.m code and save the result as Cmax.mat 
%%%%%%% parameter set as follows: deltas = 0.01; k=1;u=0;
%etaVec=0.2:1:10; lambda=2; beta=1;alpha = 1; gamma=0.62;T=10;

load('Cmax.mat')

figure
%% plot two examples of time course curve of progress 
subplot(1,4,1)
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

plot(1:T,WithoptActSeqMatrix(:,1),'o-','Color',greencolor(2,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,3),'o-','Color',greencolor(4,:));
hold on
plot(1:T,WithoptActSeqMatrix(:,7),'o-','Color',greencolor(6,:));
hold off

legend('C_{max}=0.2','C_{max}=2.2','C_{max}=6.2','Interpreter','tex')
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
legend boxoff
xlim([1,T])
xticks(1:T)
box off
set(gca,'TickDir','out');

%% days of procrastination over Cmax  
subplot(1,4,2)
numberzero = nan(1,length(etaVec));
for Idx = 1:length(etaVec)
    numberzero(Idx)= nnz(~WithoptActSeqMatrix(:,Idx));
end
plot(etaVec,numberzero,'ko-')
xlabel('maximum cost C_{max}','Interpreter','tex')
ylabel('days of procrastination')
ylim([0 T])
set(gca, 'XTick', etaVec(1:2:end));
%set(gca, 'yTick', 0:1:T);
xlim([etaVec(1),etaVec(end)])
box off
set(gca,'TickDir','out');

%% final proportional completed over Cmax  
subplot(1,4,3)
plot(etaVec,U_task_with,'ko-');
xlabel('maximum cost C_{max}','Interpreter','tex')
ylabel('final proportion completed')
ylim([0 1])
set(gca, 'XTick', etaVec(1:2:end));
xlim([etaVec(1),etaVec(end)])
box off
set(gca,'TickDir','out');

%% plot the total cost
subplot(1,4,4)
threecolor = [27,158,119;217,95,2;117,112,179]/255;
plot(etaVec,U_task_with,'o-','Color',threecolor(1,:))
hold on
plot(etaVec,cost_with,'o-','Color',threecolor(2,:));
hold on
plot(etaVec,Utotal_with,'o-','Color',threecolor(3,:))
xlabel('maximum cost C_{max}','Interpreter','tex')
ylabel('utility')
legend('final reward','total cost','net utility')
legend boxoff
ylim([0 1])
set(gca, 'XTick', etaVec(1:2:end));
xlim([etaVec(1),etaVec(end)])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])

% %% plot the U_net, R(s_T+1), total cost
% subplot(1,4,4)
% bluecolor = [158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
% % %greencolor = [153,216,201;102,194,164;65,174,118;35,139,69;0,109,44;0,68,27]/255;
% % greycolor =[189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
% % greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 
% % 
% % plot(lambdaVec,U_task_with,'o-','Color',bluecolor(3,:));
% % hold on 
% % plot(lambdaVec,cost_with,'o-','Color',redcolor(3,:));
% % hold on
% plot(etaVec,Utotal_with,'ko-');
% 
% %legend('utility of the task reward','total cost','net utility')
% xlabel('maximum cost C_{max}','Interpreter','tex')
% %ylabel('utility')
% ylabel('net utility')
% ylim([0 1])
% xlim([etaVec(1),etaVec(end)])
% set(gca, 'XTick', etaVec(1:2:end));
