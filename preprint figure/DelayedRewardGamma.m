deltas = 0.01; k=1;u=0;
eta=7; lambda=3;
alpha=1; beta=1;
deltagamma=0.1;
gammaVec = 0:deltagamma:1;
T = 10;


WithoptActSeqMatrix=nan(T,length(gammaVec));
U_task_with=nan(1,length(gammaVec));
cost_with=nan(1,length(gammaVec));
Utotal_with=nan(1,length(gammaVec));
meanWithoptActSeqMatrix = nan(1,length(gammaVec));

for gammaIdx=1:length(gammaVec)
    [WithoptActSeqMatrix(:,gammaIdx), ]=OptActStateSeq( deltas,T,k,u,eta,lambda,gammaVec(gammaIdx),alpha,beta);    
end

for gammaIdx=1:length(gammaVec)  
    meanWithoptActSeqMatrix(gammaIdx) = nanmean(WithoptActSeqMatrix(:,gammaIdx));
    U_task_with(gammaIdx)=alpha*(k*nansum(WithoptActSeqMatrix(:,gammaIdx)))^beta;
    cost_with(gammaIdx)=nansum(eta*WithoptActSeqMatrix(:,gammaIdx).^lambda);
    Utotal_with(gammaIdx) = U_task_with(gammaIdx)-cost_with(gammaIdx);
end

save('gamma.mat')
% %%%%%%% plot effort as a function of time in different discount rate conditions
% figure(1)
% bluecolor = [158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
% %greencolor = [153,216,201;102,194,164;65,174,118;35,139,69;0,109,44;0,68,27]/255;
% greycolor =[189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
% greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 
% 
% for gammaIdx=1:length(gammaVec)
%     plot(1:T,WithoptActSeqMatrix(:,gammaIdx),'o-','Color',greencolor(gammaIdx,:));
%     hold on
%     legendInfo{gammaIdx} = ['\gamma=',num2str(gammaVec(gammaIdx))];
% end
% 
% legend(legendInfo)
% xlabel('day #')
% ylabel('\Delta progress')
% xlim([1,T])
% set(gca, 'XTick', 1:T);
% set(gca, 'FontSize',28)
% 
% %%%%%% plot the proprotion completed
% figure
% plot(gammaVec,sum(WithoptActSeqMatrix,1),'ko-')
% xlabel('discount rate \gamma')
% ylabel('proportion completed')
% ylim([0,1])
% set(gca, 'FontSize',28)
% 
% 
% %%%%%% plot the U_net, R(s_T+1), total cost
% figure(2)
% plot(gammaVec,U_task_with,'o-','Color',bluecolor(3,:));
% hold on 
% plot(gammaVec,cost_with,'o-','Color',redcolor(3,:));
% hold on
% plot(gammaVec,Utotal_with,'o-','Color',greycolor(3,:));
% 
% legend('reward','total cost','net utility')
% xlabel('discount rate \gamma')
% ylabel('utility')
% legend boxoff
% set(gca, 'FontSize',28)
% 
% %%%%%%% days of procrastination over gamma  
% figure(3)
% numberzero = nan(1,length(gammaVec));
% for gammaIdx = 1:length(gammaVec)
%     numberzero(gammaIdx)= nnz(~WithoptActSeqMatrix(:,gammaIdx));
% end
% plot(gammaVec,numberzero,'ko-')
% xlabel('discount rate \gamma')
% ylabel('days of procrastination')
% ylim([0 T])
% set(gca,'FontSize',28);
% 
% %%%%%%%% plot average effort over gamma 
% % figure
% % plot(gammaVec,meanWithoptActSeqMatrix,'bo-')
% % xlabel('\gamma')
% % ylabel('average effort per day')
