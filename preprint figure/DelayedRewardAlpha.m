deltas = 0.01; k=1;u=0;
eta=10; lambda=3; beta=2;
gamma=0.4;
deltaalpha = 0.1;
alphaVec = 0:deltaalpha:1;
T = 10;

WithoptActSeqMatrix=nan(T,length(alphaVec));
U_task_with=nan(1,length(alphaVec));
cost_with=nan(1,length(alphaVec));
Utotal_with=nan(1,length(alphaVec));
finalprop = nan(1,length(alphaVec));

for alphaIdx=1:length(alphaVec)
    [WithoptActSeqMatrix(:,alphaIdx), ]=OptActStateSeq( deltas,T,k,u,eta,lambda,gamma,alphaVec(alphaIdx),beta);  
    U_task_with(alphaIdx)=alphaVec(alphaIdx)*(k*nansum(WithoptActSeqMatrix(:,alphaIdx)))^beta;
    cost_with(alphaIdx)=nansum(eta*WithoptActSeqMatrix(:,alphaIdx).^lambda);
    Utotal_with(alphaIdx) = U_task_with(alphaIdx)-cost_with(alphaIdx);
    finalprop(alphaIdx) = nansum(WithoptActSeqMatrix(:,alphaIdx));
end

save('alpha.mat')
  



%%%%%%% plot effort as a function of time in different discount rate conditions
% figure(1)
% % bluecolor = [158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
% % %greencolor = [153,216,201;102,194,164;65,174,118;35,139,69;0,109,44;0,68,27]/255;
% % greycolor =[189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
% % greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% % redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 
% 
% for alphaIdx=1:length(alphaVec)
%     plot(1:T,WithoptActSeqMatrix(:,alphaIdx),'o-','Color',greencolor(alphaIdx,:));
%     hold on
%     legendInfo{alphaIdx} = ['\alpha=',num2str(alphaVec(alphaIdx))];
% end
% 
% legend(legendInfo)
% xlabel('time t')
% ylabel('effort a')
% xlim([1,T])
% set(gca, 'XTick', 1:T);
% set(gca, 'FontSize',28)

%%%%%% plot the U_net, R(s_T+1), total cost
figure(2)
plot(alphaVec,U_task_with,'o-','Color',bluecolor(3,:));
hold on 
plot(alphaVec,cost_with,'o-','Color',redcolor(3,:));
hold on
plot(alphaVec,Utotal_with,'o-','Color',greycolor(3,:));

legend('reward','total cost','net utility')
xlabel('reward magnitude \alpha')
ylabel('utility')
legend boxoff
set(gca, 'FontSize',28)

%%%%%%% days of procrastination over alpha  
figure(3)
numberzero = nan(1,length(alphaVec));
for alphaIdx = 1:length(alphaVec)
    numberzero(alphaIdx)= nnz(~WithoptActSeqMatrix(:,alphaIdx));
end
plot(alphaVec,numberzero,'ko-')
xlabel('reward magnitude \alpha')
ylabel('days of procrastination')
ylim([0 T])
set(gca,'FontSize',28);

% %% when beta>lambda, bifurcation, average effort over alpha
% deltas = 0.01; k=1;u=0;
% eta=5; lambda=3; betaVec=4;
% alphaVec= 0.1:0.1:1; 
% gamma=1;
% T=4;
% WithoptActSeqMatrix=nan(T,length(alphaVec),length(betaVec));
% meanActSeqMatrix = nan(length(betaVec),length(alphaVec));
% U_task_with=nan(length(betaVec),length(alphaVec));
% cost_with=nan(length(betaVec),length(alphaVec));
% Utotal_with=nan(length(betaVec),length(alphaVec));
% 
% for betaIdx = 1:length(betaVec)
% for alphaIdx = 1:length(alphaVec)
%     [WithoptActSeqMatrix(:,alphaIdx,betaIdx), ]=OptActStateSeq( deltas,T,k,u,eta,lambda,gamma,alphaVec(alphaIdx),betaVec(betaIdx)); 
%     U_task_with(betaIdx,alphaIdx)=alphaVec(alphaIdx)*(nansum(WithoptActSeqMatrix(:,alphaIdx,betaIdx)))^betaVec(betaIdx);
%     cost_with(betaIdx,alphaIdx)=nansum(eta*WithoptActSeqMatrix(:,alphaIdx,betaIdx).^lambda);
%     Utotal_with(betaIdx,alphaIdx) = U_task_with(betaIdx,alphaIdx)-cost_with(betaIdx,alphaIdx);
%     meanActSeqMatrix(betaIdx,alphaIdx) = nanmean(WithoptActSeqMatrix(:,alphaIdx,betaIdx));
% end
% end
%  
% figure
% plot(alphaVec,meanActSeqMatrix,'-o','Color',greencolor(end,:))
% xlabel('\alpha')
% xlim([min(alphaVec),max(alphaVec)])
% ylabel('average effort per day')
% set(gca, 'FontSize',28,'xTick',alphaVec(1:2:end))
% 
% figure
% plot(alphaVec,k*T*meanActSeqMatrix,'-o','Color',greencolor(end,:))
% xlabel('\alpha')
% xlim([min(alphaVec),max(alphaVec)])
% ylabel('final proportion completed')
% set(gca, 'FontSize',28,'xTick',alphaVec(1:2:end))
% 
% figure
% plot(alphaVec,U_task_with,'-o','Color',bluecolor(3,:))
% hold on
% plot(alphaVec,cost_with,'-o','Color',redcolor(3,:))
% hold on
% plot(alphaVec,Utotal_with,'-o','Color',greycolor(3,:))
% legend('reward','total cost','net utility')
% xlabel('\alpha')
% ylabel('utility')
% xlim([min(alphaVec),max(alphaVec)])
% set(gca, 'FontSize',28,'xTick',alphaVec(1:2:end))
% 
% %% %%%%%%%%%%%%%%% alpha_threshould for different gammas 
% deltas = 0.01; k=1;u=0;
% eta=5; lambda=3; betaVec=4;
% alphaVec= 0.5:0.1:1; 
% %gammaVec = 0.6:0.2:1;
% gammaVec = 0.8;
% T=4;
% 
% WithoptActSeqMatrix=nan(T,length(alphaVec),length(gammaVec));
% meanActSeqMatrix = nan(length(gammaVec),length(alphaVec));
% U_task_with=nan(length(gammaVec),length(alphaVec));
% cost_with=nan(length(gammaVec),length(alphaVec));
% Utotal_with=nan(length(gammaVec),length(alphaVec));
% 
% for gammaIdx = 1:length(gammaVec)
% for alphaIdx=1:length(alphaVec)
%     [WithoptActSeqMatrix(:,alphaIdx,gammaIdx), ]=OptActStateSeq( deltas,T,k,u,eta,lambda,gammaVec(gammaIdx),alphaVec(alphaIdx),beta); 
%     U_task_with(gammaIdx,alphaIdx)=alphaVec(alphaIdx)*(k*nansum(WithoptActSeqMatrix(:,alphaIdx,gammaIdx)))^beta;
%     cost_with(gammaIdx,alphaIdx)=nansum(eta*WithoptActSeqMatrix(:,alphaIdx,gammaIdx).^lambda);
%     Utotal_with(gammaIdx,alphaIdx) = U_task_with(gammaIdx,alphaIdx)-cost_with(gammaIdx,alphaIdx);
%     meanActSeqMatrix(gammaIdx,alphaIdx) = nanmean(WithoptActSeqMatrix(:,alphaIdx,gammaIdx));
% end
% end
% 
% figure
% for gammaIdx = 1:length(gammaVec)
%     plot(alphaVec,meanActSeqMatrix(gammaIdx,:),'-o','Color',bluecolor(2*gammaIdx,:))
%     hold on 
%     legendInfo{gammaIdx} = ['\gamma=',num2str(gammaVec(gammaIdx))];
% end
% legend(legendInfo)
% xlabel('\alpha')
% xlim([min(alphaVec),max(alphaVec)])
% ylabel('average effort per day')
% set(gca, 'FontSize',28,'xTick',alphaVec(1:2:end))
% 
% figure
% for gammaIdx = 1:length(gammaVec)
%     plot(alphaVec,k*T*meanActSeqMatrix(gammaIdx,:),'-o','Color',greencolor(2*gammaIdx,:))
%     hold on 
%     legendInfo{gammaIdx} = ['\gamma=',num2str(gammaVec(gammaIdx))];
% end
% legend(legendInfo)
% ylabel('final proportion completed')
% xlim([min(alphaVec),max(alphaVec)])
% xlabel('\alpha')
% set(gca, 'FontSize',28,'xTick',alphaVec(1:2:end))
% 
% figure
% for gammaIdx =1:length(gammaVec)
% plot(alphaVec,U_task_with(gammaIdx,:),'o-','Color',bluecolor(gammaIdx,:));
% hold on 
% plot(alphaVec,cost_with(gammaIdx,:),'o-','Color',redcolor(gammaIdx,:));
% hold on
% plot(alphaVec,Utotal_with(gammaIdx,:),'o-','Color',greycolor(gammaIdx,:));
% end
% xlim([min(alphaVec),max(alphaVec)])
% xlabel('\alpha')
% ylabel('utility')
% set(gca, 'FontSize',28,'xTick',alphaVec(1:2:end))
% legend('reward','total cost','net utility')
% 

% figure
% for betaIdx = 1:length(betaVec)
%     plot(alphaVec,meanActSeqMatrix(betaIdx,:),'-o','Color',bluecolor(betaIdx+3,:))
%     hold on 
%     legendInfo{betaIdx} = ['\beta=',num2str(betaVec(betaIdx))];
% end
% legend(legendInfo)
% xlabel('\alpha')
% xlim([min(alphaVec),max(alphaVec)])
% ylabel('average effort per day')
% set(gca, 'FontSize',28,'xTick',alphaVec)

% figure
% subplot(1,3,1)
% for betaIdx = 1:length(betaVec)
% plot(alphaVec,U_task_with(betaIdx,:),'-o','Color',greencolor(betaIdx+3,:))
% hold on
% end
% xlabel('\alpha')
% ylabel('reward')
% xlim([min(alphaVec),max(alphaVec)])
% set(gca, 'FontSize',28,'xTick',alphaVec)
% subplot(1,3,2)
% for betaIdx = 1:length(betaVec)
% plot(alphaVec,cost_with(betaIdx,:),'-o','Color',redcolor(betaIdx+3,:))
% hold on
% end
% xlabel('\alpha')
% ylabel('total cost')
% xlim([min(alphaVec),max(alphaVec)])
% set(gca, 'FontSize',28,'xTick',alphaVec)
% subplot(1,3,3)
% for betaIdx = 1:length(betaVec)
% plot(alphaVec,Utotal_with(betaIdx,:),'-o','Color',greycolor(betaIdx+3,:))
% hold on
% legendInfo{betaIdx} = ['\beta=',num2str(betaVec(betaIdx))];
% end
% legend(legendInfo)
% xlabel('\alpha')
% ylabel('net utility')
% xlim([min(alphaVec),max(alphaVec)])
% set(gca, 'FontSize',28,'xTick',alphaVec)

% %% change gammma here
% deltas = 0.01; k=1;u=0;
% eta=0.4; lambda=3; beta=4;
% alphaVec= 0.1:0.2:1; 
% gammaVec = 0:0.3:1; Tvec =1:10; 
% % when gamma=1
% WithoptActSeqMatrix=nan(Tvec(end),Tvec(end),length(gammaVec));
% meanActSeqMatrix = nan(length(gammaVec),length(Tvec));
% U_task_with=nan(length(gammaVec),length(Tvec));
% cost_with=nan(length(gammaVec),length(Tvec));
% Utotal_with=nan(length(gammaVec),length(Tvec));
% 
% for gammaIdx = 1:length(gammaVec)
% for TvecIdx=1:length(Tvec)
%     [WithoptActSeqMatrix(1:Tvec(TvecIdx),TvecIdx,gammaIdx), ]=OptActStateSeq( deltas,Tvec(TvecIdx),k,u,eta,lambda,gammaVec(gammaIdx),alphaVec(end),beta); 
%     U_task_with(gammaIdx,TvecIdx)=alphaVec(end)*(nansum(WithoptActSeqMatrix(:,TvecIdx,gammaIdx)))^beta;
%     cost_with(gammaIdx,TvecIdx)=nansum(eta*WithoptActSeqMatrix(:,TvecIdx,gammaIdx).^lambda);
%     Utotal_with(gammaIdx,TvecIdx) = U_task_with(gammaIdx,TvecIdx)-cost_with(gammaIdx,TvecIdx);
%     meanActSeqMatrix(gammaIdx,TvecIdx) = nanmean(WithoptActSeqMatrix(:,TvecIdx,gammaIdx));
% end
% end
% 
% figure
% for Tidx =1:length(Tvec)-1
%     plot(1:Tidx,WithoptActSeqMatrix(1:Tidx,Tidx,2),'-o','Color',bluecolor(Tidx,:))
%     hold on
% end
% figure
% for gammaIdx = 1:length(gammaVec)
%     plot(Tvec,meanActSeqMatrix(gammaIdx,:),'-o','Color',bluecolor(gammaIdx+3,:))
%     hold on 
%     legendInfo{gammaIdx} = ['\gamma=',num2str(gammaVec(gammaIdx))];
% end
% legend(legendInfo)
% xlabel('total time T')
% ylabel('average effort per day')
% xlim([Tvec(1),Tvec(end)])
% set(gca, 'XTick', Tvec);
% set(gca, 'FontSize',28)
% 
% figure
% subplot(1,3,1)
% for gammaIdx = 1:length(gammaVec)
% 
%     plot(Tvec,U_task_with(gammaIdx,:),'-o','Color',greencolor(gammaIdx+3,:))
%     hold on     
%     legendInfo{gammaIdx} = ['\gamma=',num2str(gammaVec(gammaIdx))];
% end
% legend(legendInfo)
% xlabel('total time T')
% ylabel('reward')
% xlim([Tvec(1),Tvec(end)])
% set(gca, 'XTick', Tvec);
% set(gca, 'FontSize',28)
% subplot(1,3,2)
% for gammaIdx = 1:length(gammaVec)
% 
%     plot(Tvec,cost_with(gammaIdx,:),'-o','Color',redcolor(gammaIdx+3,:))
%     hold on     
%     legendInfo{gammaIdx} = ['\gamma=',num2str(gammaVec(gammaIdx))];
% end
% legend(legendInfo)
% xlabel('total time T')
% ylabel('total cost')
% xlim([Tvec(1),Tvec(end)])
% set(gca, 'XTick', Tvec);
% set(gca, 'FontSize',28)
% subplot(1,3,3)
% for gammaIdx = 1:length(gammaVec)
% 
%     plot(Tvec,Utotal_with(gammaIdx,:),'-o','Color',greycolor(gammaIdx+3,:))
%     hold on     
%     legendInfo{gammaIdx} = ['\gamma=',num2str(gammaVec(gammaIdx))];
% end
% legend(legendInfo)
% xlabel('total time T')
% ylabel('net utility')
% xlim([Tvec(1),Tvec(end)])
% set(gca, 'XTick', Tvec);
% set(gca, 'FontSize',28)

% %% T=6, gamma=0.3, for  alpha<c1 and alpha>c1, the average effort per day is 1/T. Let?s see the effort over t for different , and its net utility function. 
% deltas = 0.01; k=1;u=0;
% eta=0.4; lambda=3; beta=0.5;
% alphaVec= 0.1:0.2:0.9; 
% gammaVec = 0:0.2:1; T=6;% gamma = 0.3;
% 
% WithoptActSeqMatrix=nan(length(alphaVec),T);
% U_task_with=nan(1,length(alphaVec));
% cost_with=nan(1,length(alphaVec));
% Utotal_with=nan(1,length(alphaVec));
% 
% %% change alpha here
% for alphaIdx = 1:length(alphaVec)
%     [WithoptActSeqMatrix(alphaIdx,:), ]=OptActStateSeq( deltas,T,k,u,eta,lambda,gamma,alphaVec(alphaIdx),beta); 
%     U_task_with(alphaIdx)=alphaVec(alphaIdx)*(nansum(WithoptActSeqMatrix(alphaIdx,:)))^beta;
%     cost_with(alphaIdx)=nansum(eta*WithoptActSeqMatrix(alphaIdx,:).^lambda);
%     Utotal_with(alphaIdx) = U_task_with(alphaIdx)-cost_with(alphaIdx);
% end
% 
% figure
% for alphaIdx = 1:length(alphaVec)
%     plot(1:T,WithoptActSeqMatrix(alphaIdx,:),'-o','Color',bluecolor(alphaIdx+3,:))
%     hold on
% end
% xlabel('time t')
% ylabel('effort a')
% 
% figure
% plot(alphaVec,U_task_with,'-o','Color',greencolor(end-2,:))
% hold on
% plot(alphaVec,cost_with,'-o','Color',redcolor(end-2,:))
% hold on
% plot(alphaVec,Utotal_with,'-o','Color',greycolor(end-2,:))
% legend('reward','total cost','net utility')
% xlabel('\alpha')
% ylabel('utility')