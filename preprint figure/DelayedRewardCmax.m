%%%% orignal file name DelayedRewardEta.m in folder FirstDraftPrograms

%%%%%% change over c_1 (eta)
deltas = 0.01; k=1;u=0;
etaVec=0.2:1:10; lambda=2; beta=1;
alpha = 1; 
gamma=0.62;
T=10;
WithoptActSeqMatrix=nan(T,length(etaVec));
meanActSeqMatrix = nan(1,length(etaVec));
U_task_with=nan(1,length(etaVec));
cost_with=nan(1,length(etaVec));
Utotal_with=nan(1,length(etaVec));

for etaIdx = 1:length(etaVec)
    [WithoptActSeqMatrix(:,etaIdx), ]=OptActStateSeq( deltas,T,k,u,etaVec(etaIdx),lambda,gamma,alpha,beta); 
    U_task_with(etaIdx)=alpha*(nansum(WithoptActSeqMatrix(:,etaIdx)))^beta;
    cost_with(etaIdx)=nansum(etaVec(etaIdx)*WithoptActSeqMatrix(:,etaIdx).^lambda);
    Utotal_with(etaIdx) = U_task_with(etaIdx)-cost_with(etaIdx);
    meanActSeqMatrix(etaIdx) = nanmean(WithoptActSeqMatrix(:,etaIdx));
end

save('Cmax.mat')
% %%%%%%% plot effort as a function of time in different discount rate conditions
% figure(1)
% bluecolor = [158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
% greycolor =[189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
% greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 
% 
% for etaIdx=1:length(etaVec)
%     plot(1:T,WithoptActSeqMatrix(:,etaIdx),'o-','Color',greencolor(etaIdx,:));
%     hold on
%     legendInfo{etaIdx} = ['c_1=',num2str(etaVec(etaIdx))];
% end
% 
% legend(legendInfo)
% xlabel('time t')
% ylabel('effort a')
% xlim([1,T])
% set(gca, 'XTick', 1:T);
% set(gca, 'FontSize',28)
% 
% %%%%%% plot the U_net, R(s_T+1), total cost
% figure(2)
% plot(etaVec,U_task_with,'o-','Color',bluecolor(3,:));
% hold on 
% plot(etaVec,cost_with,'o-','Color',redcolor(3,:));
% hold on
% plot(etaVec,Utotal_with,'o-','Color',greycolor(3,:));
% 
% legend('reward','total cost','net utility')
% xlabel('cost magnitude c_1')
% ylabel('utility')
% legend boxoff
% set(gca, 'FontSize',28)

%%%%%%% days of procrastination over alpha  
% figure(3)
% numberzero = nan(1,length(etaVec));
% for etaIdx = 1:length(etaVec)
%     numberzero(etaIdx)= nnz(~WithoptActSeqMatrix(:,etaIdx));
% end
% plot(etaVec,numberzero,'ko-')
% xlabel('cost magnitude c_1')
% ylabel('days of procrastination')
% ylim([0 T])
% set(gca,'FontSize',28);

% greycolor =[255,255,255;240,240,240;217,217,217;189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
% greencolor = [247,252,245;229,245,224;199,233,192;161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% redcolor = [255,245,240;254,224,210;252,187,161;252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 
% bluecolor = [247,251,255;247,251,255;222,235,247;198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;

% figure
% plot(etaVec,meanActSeqMatrix,'-o','Color',bluecolor(end-1,:))
% xlabel('c_1')
% xlim([min(etaVec),max(etaVec)])
% ylabel('average effort per day')
% set(gca, 'FontSize',28,'xTick',etaVec)
% 
% figure
% plot(etaVec,U_task_with,'-o','Color',greencolor(end-1,:))
% hold on
% plot(etaVec,cost_with,'-o','Color',redcolor(end-1,:))
% hold on
% plot(etaVec,Utotal_with,'-o','Color',greycolor(end-1,:))
% legend('reward','total cost','net utility')
% xlabel('c_1')
% ylabel('utility')
% xlim([min(etaVec),max(etaVec)])
% set(gca, 'FontSize',28,'xTick',etaVec)
% 
% %% c_1_threshould for different gammas 
% deltas = 0.01; k=1;u=0;
% etaVec=1:10; lambda=3; beta=4;
% alpha = 0.4; 
% gammaVec=0.2:0.2:1;
% T=4;
% 
% WithoptActSeqMatrix=nan(T,length(etaVec),length(gammaVec));
% meanActSeqMatrix = nan(length(etaVec),length(kVec));
% U_task_with=nan(length(gammaVec),length(etaVec));
% cost_with=nan(length(gammaVec),length(etaVec));
% Utotal_with=nan(length(gammaVec),length(etaVec));
% 
% for gammaIdx = 1:length(gammaVec)
% for etaIdx=1:length(etaVec)
%     [WithoptActSeqMatrix(:,etaIdx,gammaIdx), ]=OptActStateSeq( deltas,T,k,u,etaVec(etaIdx),lambda,gammaVec(gammaIdx),alpha,beta); 
%     U_task_with(gammaIdx,etaIdx)=alpha*(k*nansum(WithoptActSeqMatrix(:,etaIdx,gammaIdx)))^beta;
%     cost_with(gammaIdx,etaIdx)=nansum(etaVec(etaIdx)*WithoptActSeqMatrix(:,etaIdx,gammaIdx).^lambda);
%     Utotal_with(gammaIdx,etaIdx) = U_task_with(gammaIdx,etaIdx)-cost_with(gammaIdx,etaIdx);
%     meanActSeqMatrix(gammaIdx,etaIdx) = nanmean(WithoptActSeqMatrix(:,etaIdx,gammaIdx));
% end
% end
% 
% figure
% for gammaIdx = 1:length(gammaVec)
%     plot(etaVec,meanActSeqMatrix(gammaIdx,:),'-o','Color',bluecolor(gammaIdx+3,:))
%     hold on 
%     legendInfo{gammaIdx} = ['\gamma=',num2str(gammaVec(gammaIdx))];
% end
% legend(legendInfo)
% xlabel('c_1')
% xlim([min(etaVec),max(etaVec)])
% ylabel('average effort per day')
% set(gca, 'FontSize',28,'xTick',etaVec)
% 
% figure
% for gammaIdx =1:length(gammaVec)
% plot(etaVec,U_task_with(gammaIdx,:),'o-','Color',greencolor(gammaIdx+3,:));
% hold on 
% plot(etaVec,cost_with(gammaIdx,:),'o-','Color',redcolor(gammaIdx+3,:));
% hold on
% plot(etaVec,Utotal_with(gammaIdx,:),'o-','Color',greycolor(gammaIdx+3,:));
% end
% xlim([min(etaVec),max(etaVec)])
% xlabel('c_1')
% ylabel('utility')
% set(gca, 'FontSize',28,'xTick',etaVec)
% legend('reward','total cost','net utility')