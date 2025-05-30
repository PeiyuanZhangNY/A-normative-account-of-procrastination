close all;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

deltas = 0.01; k=1;u=0;
eta=7; lambda=3;
alpha=1; 
plotgammaVec = [0.1,0.4,0.35];
%gamma = 0.35;
%gamma = 0.4; 
%gamma = 0.1;
T = 10; 
betaVec = 0.5:0.5:3.5; 
%betaVec = 1:0.5:4;
plotbetaIndex = [1,4,7];

WithoptActSeqMatrix=nan(length(plotgammaVec),length(betaVec),T);
finalprop = nan(length(plotgammaVec),length(betaVec));
U_task_with=nan(length(plotgammaVec),length(betaVec));
cost_with=nan(length(plotgammaVec),length(betaVec));
Utotal_with=nan(length(plotgammaVec),length(betaVec));

for gammaIdx = 1:length(plotgammaVec)
for betaIdx=1:length(betaVec)
    
    [WithoptActSeqMatrix(gammaIdx,betaIdx,:), ]=OptActStateSeq( deltas,T,k,u,eta,lambda,plotgammaVec(gammaIdx),alpha,betaVec(betaIdx)); 
    temp = reshape(WithoptActSeqMatrix(gammaIdx,betaIdx,:),[1,T]);
    U_task_with(gammaIdx,betaIdx)=alpha*(nansum(temp)^betaVec(betaIdx));
    cost_with(gammaIdx,betaIdx)=nansum(eta*temp.^lambda);
    Utotal_with(gammaIdx,betaIdx) = U_task_with(gammaIdx,betaIdx)-cost_with(gammaIdx,betaIdx);
    finalprop(gammaIdx,betaIdx) = nansum(temp);

end
end

bluecolor = [158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
greycolor =[189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 

%%%%%% for gamma= 0.1
%%%%% plot effort as a function of time in different beta conditions
figure
subplot(1,4,1)
plot(1:T,reshape(WithoptActSeqMatrix(1,plotbetaIndex(1),:),[1,T]),'o-','Color',greencolor(2,:));
hold on
plot(1:T,reshape(WithoptActSeqMatrix(1,plotbetaIndex(2),:),[1,T]),'o-','Color',greencolor(4,:));
hold on
plot(1:T,reshape(WithoptActSeqMatrix(1,plotbetaIndex(3),:),[1,T]),'o-','Color',greencolor(6,:));
hold off

legend('\beta=0.5','\beta=2','\beta=3.5')
legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
set(gca, 'XTick', 1:T);
box off
set(gca,'TickDir','out');

%%%%%%% days of procrastination over beta
subplot(1,4,2)
numberzero = nan(1,length(betaVec));
for betaIdx = 1:length(betaVec)
    numberzero(betaIdx)= nnz(~WithoptActSeqMatrix(1,betaIdx,:));
end
plot(betaVec,numberzero,'ko-')
xlabel('level of perfectionism \beta')
ylabel('days of procrastination')
ylim([0 T])
xlim([betaVec(1),betaVec(end)])
set(gca, 'XTick', betaVec(1):1:betaVec(end));
box off
set(gca,'TickDir','out');

%%%%%%% final proportion completed
subplot(1,4,3)
plot(betaVec,finalprop(1,:),'ko-')
ylim([0 1])
xlim([betaVec(1),betaVec(end)])
xlabel('level of perfectionism \beta')
ylabel('final proportion completed')
set(gca, 'XTick', betaVec(1):1:betaVec(end));
box off
set(gca,'TickDir','out');

% %%%%%% plot total cost
subplot(1,4,4)
threecolor = [27,158,119;217,95,2;117,112,179]/255;
plot(betaVec,U_task_with(1,:),'o-','Color',threecolor(1,:));
hold on 
plot(betaVec,cost_with(1,:),'o-','Color',threecolor(2,:));
hold on
plot(betaVec,Utotal_with(1,:),'o-','Color',threecolor(3,:));
legend('final reward','total cost','net utility')
legend boxoff

xlabel('level of perfectionism \beta')
ylabel('utility')
set(gca, 'XTick', betaVec(1):1:betaVec(end));
ylim([0 1])
xlim([betaVec(1),betaVec(end)])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])

%%%%%% for gamma= 0.4
%%%%% plot effort as a function of time in different beta conditions
figure
subplot(1,4,1)
plot(1:T,reshape(WithoptActSeqMatrix(2,plotbetaIndex(1),:),[1,T]),'o-','Color',greencolor(2,:));
hold on
plot(1:T,reshape(WithoptActSeqMatrix(2,plotbetaIndex(2),:),[1,T]),'o-','Color',greencolor(4,:));
hold on
plot(1:T,reshape(WithoptActSeqMatrix(2,plotbetaIndex(3),:),[1,T]),'o-','Color',greencolor(6,:));
hold off

legend('\beta=0.5','\beta=2','\beta=3.5')
legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
set(gca, 'XTick', 1:T);
box off
set(gca,'TickDir','out');

%%%%%%% days of procrastination over beta
subplot(1,4,2)
numberzero = nan(1,length(betaVec));
for betaIdx = 1:length(betaVec)
    numberzero(betaIdx)= nnz(~WithoptActSeqMatrix(2,betaIdx,:));
end
plot(betaVec,numberzero,'ko-')
xlabel('level of perfectionism \beta')
ylabel('days of procrastination')
ylim([0 T])
xlim([betaVec(1),betaVec(end)])
set(gca, 'XTick', betaVec(1):1:betaVec(end));
box off
set(gca,'TickDir','out');

%%%%%%% final proportion completed
subplot(1,4,3)
plot(betaVec,finalprop(2,:),'ko-')
ylim([0 1])
xlim([betaVec(1),betaVec(end)])
xlabel('level of perfectionism \beta')
ylabel('final proportion completed')
set(gca, 'XTick', betaVec(1):1:betaVec(end));
box off
set(gca,'TickDir','out');

% %%%%%% plot the U_net, R(s_T+1), total cost
subplot(1,4,4)
threecolor = [27,158,119;217,95,2;117,112,179]/255;
plot(betaVec,U_task_with(2,:),'o-','Color',threecolor(1,:));
hold on 
plot(betaVec,cost_with(2,:),'o-','Color',threecolor(2,:));
hold on
plot(betaVec,Utotal_with(2,:),'o-','Color',threecolor(3,:));
legend('final reward','total cost','net utility')
legend boxoff
xlabel('level of perfectionism \beta')
ylabel('utility')
set(gca, 'XTick', betaVec(1):1:betaVec(end));
ylim([0 1])
xlim([betaVec(1),betaVec(end)])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])
%%%%%% for gamma= 0.35
%%%%% plot effort as a function of time in different beta conditions
figure
subplot(1,4,1)
plot(1:T,reshape(WithoptActSeqMatrix(3,plotbetaIndex(1),:),[1,T]),'o-','Color',greencolor(2,:));
hold on
plot(1:T,reshape(WithoptActSeqMatrix(3,plotbetaIndex(2),:),[1,T]),'o-','Color',greencolor(4,:));
hold on
plot(1:T,reshape(WithoptActSeqMatrix(3,plotbetaIndex(3),:),[1,T]),'o-','Color',greencolor(6,:));
hold off

legend('\beta=0.5','\beta=2','\beta=3.5')
legend boxoff
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
set(gca, 'XTick', 1:T);
box off
set(gca,'TickDir','out');

%%%%%%% days of procrastination over beta
subplot(1,4,2)
numberzero = nan(1,length(betaVec));
for betaIdx = 1:length(betaVec)
    numberzero(betaIdx)= nnz(~WithoptActSeqMatrix(3,betaIdx,:));
end
plot(betaVec,numberzero,'ko-')
xlabel('level of perfectionism \beta')
ylabel('days of procrastination')
ylim([0 T])
xlim([betaVec(1),betaVec(end)])
set(gca, 'XTick', betaVec(1):1:betaVec(end));
box off
set(gca,'TickDir','out');

%%%%%%% final proportion completed
subplot(1,4,3)
plot(betaVec,finalprop(3,:),'ko-')
ylim([0 1])
xlim([betaVec(1),betaVec(end)])
xlabel('level of perfectionism \beta')
ylabel('final proportion completed')
set(gca, 'XTick', betaVec(1):1:betaVec(end));
box off
set(gca,'TickDir','out');

% %%%%%% plot the U_net, R(s_T+1), total cost
subplot(1,4,4)
threecolor = [27,158,119;217,95,2;117,112,179]/255;
plot(betaVec,U_task_with(3,:),'o-','Color',threecolor(1,:));
hold on 
plot(betaVec,cost_with(3,:),'o-','Color',threecolor(2,:));
hold on
plot(betaVec,Utotal_with(3,:),'o-','Color',threecolor(3,:));
legend('final reward','total cost','net utility')
legend boxoff
xlabel('level of perfectionism \beta')
ylabel('utility')
set(gca, 'XTick', betaVec(1):1:betaVec(end));
ylim([0 1])
xlim([betaVec(1),betaVec(end)])
box off
set(gca,'TickDir','out');

x0=10;
y0=10;
width=1200;
height=250;
set(gcf,'position',[x0,y0,width,height])
%% average effort over gamma in different beta 
%betaVec = 0.5:0.5:3.5; 
%plotbetaIndex = [1,5,7];
gammaVec = 0:0.1:1;
WithoptActSeqMatrix=nan(T,length(gammaVec),length(betaVec));
meanActSeqMatrix = nan(length(betaVec),length(gammaVec));
U_task_with=nan(length(betaVec),length(gammaVec));
cost_with=nan(length(betaVec),length(gammaVec));
Utotal_with=nan(length(betaVec),length(gammaVec));

for betaIdx = 1:length(betaVec)
for gammaIdx=1:length(gammaVec)
    [WithoptActSeqMatrix(:,gammaIdx,betaIdx), ]=OptActStateSeq( deltas,T,k,u,eta,lambda,gammaVec(gammaIdx),alpha,betaVec(betaIdx)); 
    U_task_with(betaIdx,gammaIdx)=alpha*(k*nansum(WithoptActSeqMatrix(:,gammaIdx,betaIdx)))^betaVec(betaIdx);
    cost_with(betaIdx,gammaIdx)=nansum(eta*WithoptActSeqMatrix(:,gammaIdx,betaIdx).^lambda);
    Utotal_with(betaIdx,gammaIdx) = U_task_with(betaIdx,gammaIdx)-cost_with(betaIdx,gammaIdx);
    meanActSeqMatrix(betaIdx,gammaIdx) = nanmean(WithoptActSeqMatrix(:,gammaIdx,betaIdx));
end
end

figure
plot(gammaVec,k*T*meanActSeqMatrix(plotbetaIndex(1),:),'-o','Color',greencolor(2,:))
hold on 
plot(gammaVec,k*T*meanActSeqMatrix(plotbetaIndex(2),:),'-o','Color',greencolor(4,:))
hold on
plot(gammaVec,k*T*meanActSeqMatrix(plotbetaIndex(3),:),'-o','Color',greencolor(6,:))
hold off

  
legend('\beta=0.5','\beta=2','\beta=3.5')
legend boxoff
xlabel('discount rate \gamma')
ylabel('final proportion completed')
xlim([gammaVec(1),gammaVec(end)])
set(gca, 'XTick', gammaVec(1:2:end));
axis square
box off
set(gca,'TickDir','out');
% figure
% for betaIdx = 1:length(betaVec)
%     plot([gammaVec(3),gammaVec(5)],[meanActSeqMatrix(betaIdx,3),meanActSeqMatrix(betaIdx,5)],'o','Color',greencolor(2*betaIdx,:))
%     hold on 
%     legendInfo{betaIdx} = ['\beta=',num2str(betaVec(betaIdx))];
% end
% legend(legendInfo)
% xlabel('\gamma')
% ylabel('average effort per day')
% xlim([gammaVec(1),gammaVec(end)])
% set(gca, 'XTick', gammaVec(1:2:end));
% set(gca, 'FontSize',28)
% 
% figure
% for betaIdx = 1:length(betaVec)
%     plot(gammaVec,meanActSeqMatrix(betaIdx,:),'-o','Color',greencolor(2*betaIdx,:))
%     hold on 
%     legendInfo{betaIdx} = ['\beta=',num2str(betaVec(betaIdx))];
% end
% legend(legendInfo)
% xlabel('\gamma')
% ylabel('average effort per day')
% xlim([gammaVec(1),gammaVec(end)])
% set(gca, 'XTick', gammaVec(1:2:end));
% set(gca, 'FontSize',28)
% 

% set(gca, 'FontSize',28)
% 
% figure
% plot(gammaVec,U_task_with(end,:),'o-','Color',bluecolor(3,:));
% hold on 
% plot(gammaVec,cost_with(end,:),'o-','Color',redcolor(3,:));
% hold on
% plot(gammaVec,Utotal_with(end,:),'o-','Color',greycolor(3,:));
% legend('reward','total cost','net utility')
% xlabel('\gamma')
% ylabel('utility')
% xlim([gammaVec(1),gammaVec(end)])
% set(gca, 'XTick', gammaVec(1:2:end));
% set(gca, 'FontSize',28)
% 
% 
% %% average effort over c_1 in different beta 
% etaVec = 1:2:19;
% WithoptActSeqMatrix=nan(T,length(etaVec),length(betaVec));
% meanActSeqMatrix = nan(length(betaVec),length(etaVec));
% U_task_with=nan(length(betaVec),length(etaVec));
% cost_with=nan(length(betaVec),length(etaVec));
% Utotal_with=nan(length(betaVec),length(etaVec));
% 
% for betaIdx = 1:length(betaVec)
% for etaIdx=1:length(etaVec)
%     [WithoptActSeqMatrix(:,etaIdx,betaIdx), ]=OptActStateSeq( deltas,T,k,u,etaVec(etaIdx),lambda,gamma,alpha,betaVec(betaIdx)); 
%     U_task_with(betaIdx,etaIdx)=alpha*(k*nansum(WithoptActSeqMatrix(:,etaIdx,betaIdx)))^betaVec(betaIdx);
%     cost_with(betaIdx,etaIdx)=nansum(etaVec(etaIdx)*WithoptActSeqMatrix(:,etaIdx,betaIdx).^lambda);
%     Utotal_with(betaIdx,etaIdx) = U_task_with(betaIdx,etaIdx)-cost_with(betaIdx,etaIdx);
%     meanActSeqMatrix(betaIdx,etaIdx) = nanmean(WithoptActSeqMatrix(:,etaIdx,betaIdx));
% end
% end
% 
% figure
% for betaIdx = 1:length(betaVec)
%     plot([etaVec(4),etaVec(8)],[meanActSeqMatrix(betaIdx,4),meanActSeqMatrix(betaIdx,8)],'o','Color',greencolor(2*betaIdx,:))
%     hold on 
%     legendInfo{betaIdx} = ['\beta=',num2str(betaVec(betaIdx))];
% end
% legend(legendInfo)
% xlabel('cost magnitude c_1')
% ylabel('average effort per day')
% xlim([etaVec(1),etaVec(end)])
% set(gca, 'XTick', etaVec);
% set(gca, 'FontSize',28)
% 
% figure
% for betaIdx = 1:length(betaVec)
%     plot(etaVec,meanActSeqMatrix(betaIdx,:),'-o','Color',greencolor(2*betaIdx,:))
%     hold on 
%     legendInfo{betaIdx} = ['\beta=',num2str(betaVec(betaIdx))];
% end
% legend(legendInfo)
% xlabel('cost magnitude c_1')
% ylabel('average effort per day')
% xlim([etaVec(1),etaVec(end)])
% set(gca, 'XTick', etaVec);
% 
% figure
% for betaIdx = 1:length(betaVec)
%     plot(etaVec,k*T*meanActSeqMatrix(betaIdx,:),'-o','Color',greencolor(2*betaIdx,:))
%     hold on 
%     legendInfo{betaIdx} = ['\beta=',num2str(betaVec(betaIdx))];
% end
% legend(legendInfo)
% xlabel('cost magnitude c_1')
% ylabel('final proportion completed')
% xlim([etaVec(1),etaVec(end)])
% set(gca, 'XTick', etaVec);
% set(gca, 'FontSize',28)
% 
% figure
% plot(etaVec,U_task_with(end,:),'o-','Color',bluecolor(3,:));
% hold on 
% plot(etaVec,cost_with(end,:),'o-','Color',redcolor(3,:));
% hold on
% plot(etaVec,Utotal_with(end,:),'o-','Color',greycolor(3,:));
% legend('reward','total cost','net utility')
% xlabel('cost magnitude c_1')
% ylabel('utility')
% xlim([etaVec(1),etaVec(end)])
% set(gca, 'XTick', etaVec);
% set(gca, 'FontSize',28)
% 
