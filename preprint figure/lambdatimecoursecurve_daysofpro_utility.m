close all;
clear
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

%%%%%%%% first run DelayedRewardLambda.m code 

%%%%%%% parameters for figure
%%%%%%% lambdatimecoursecurve_daysofpro_utility.pdf
% eta=0.8; gamma = 0.2;
%%%%%%% lambdaVec = 0.2:0.4:4;
load('lambda.mat')

figure
%% plot examples of time course curve of progress 
subplot(1,4,1)
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

plot(1:T,WithoptActSeqMatrix(:,2),'o-','Color',[216,179,101]/255);
hold on
plot(1:T,WithoptActSeqMatrix(:,5),'o-','Color',[90,180,172]/255);
hold on
plot(1:T,WithoptActSeqMatrix(:,8),'o-','Color',([90,180,172]+30)/255);%greencolor(4,:));
hold off

legend('\lambda=0.6','\lambda=1.8','\lambda=3')
legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
xticks(1:T)
box off
set(gca,'TickDir','out');

%% days of procrastination over lambda  
subplot(1,4,2)
numberzero = nan(1,length(lambdaVec));
for lambdaIdx = 1:length(lambdaVec)
    numberzero(lambdaIdx)= nnz(~WithoptActSeqMatrix(:,lambdaIdx));
end
fill([min(lambdaVec) 1 1 min(lambdaVec)], [0 0 7 7], [216,179,101]/255, 'EdgeColor', 'none','FaceAlpha',0.6);
hold on
plot(lambdaVec,numberzero,'ko-')
xlabel('exponent of cost function \it\lambda','Interpreter','tex')
ylabel('days of procrastination')
ylim([0 T])
set(gca, 'XTick', lambdaVec(1:2:end));
xlim([lambdaVec(1),lambdaVec(end)])
box off
set(gca,'TickDir','out');

%% final proportion completed over lambda  
subplot(1,4,3)
fill([min(lambdaVec) 1 1 min(lambdaVec)], [0 0 7 7], [216,179,101]/255, 'EdgeColor', 'none','FaceAlpha',0.6);
hold on
plot(lambdaVec,U_task_with,'ko-');
xlabel('exponent of cost function \it\lambda','Interpreter','tex')
ylabel('final proportion completed')
ylim([0 1])
set(gca, 'XTick', lambdaVec(1:2:end));
xlim([lambdaVec(1),lambdaVec(end)])
box off
set(gca,'TickDir','out');

%% plot the total cost
subplot(1,4,4)
bluecolor = [158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
% %greencolor = [153,216,201;102,194,164;65,174,118;35,139,69;0,109,44;0,68,27]/255;
% greycolor =[189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
% greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 
threecolor = [27,158,119;217,95,2;117,112,179]/255;
% 
% plot(lambdaVec,U_task_with,'o-','Color',bluecolor(3,:));
% hold on 
% plot(lambdaVec,cost_with,'o-','Color',redcolor(3,:));
% hold on
fill([min(lambdaVec) 1 1 min(lambdaVec)], [0 0 7 7], [216,179,101]/255, 'EdgeColor', 'none','FaceAlpha',0.6);
hold on
plot(lambdaVec,U_task_with,'o-','Color',threecolor(1,:))
hold on
plot(lambdaVec,cost_with,'o-','Color',threecolor(2,:))
hold on
plot(lambdaVec,Utotal_with,'o-','Color',threecolor(3,:));
legend('','final reward','total cost','net utility')
legend boxoff
%legend('utility of the task reward','total cost','net utility')
xlabel('exponent of cost function \it\lambda','Interpreter','tex')
%ylabel('utility')
ylabel('utility')
ylim([0 1])
xlim([lambdaVec(1),lambdaVec(end)])
set(gca, 'XTick', lambdaVec(1:2:end));
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])