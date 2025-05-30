% deltas = 0.01; alpha =0.8;u=0;
% eta=55; lambda=5; beta=4.5;
% gamma=0.2; deltaT = 1;
% k=1;TVec = 2:deltaT:10;

deltas = 0.01; alpha =1;u=0;
eta=10; lambda=3; beta=2.8;
gamma=0.6; deltaT = 1;
k=1;TVec = 2:deltaT:17;

WithoptActSeqMatrix=nan(TVec(end),length(TVec));
U_task_with=nan(1,length(TVec));
cost_with=nan(1,length(TVec));
Utotal_with=nan(1,length(TVec));
finalprop = nan(1,length(TVec));

for TIdx=1:length(TVec)
    T = TVec(TIdx);
    [WithoptActSeqMatrix(1:T,TIdx), ]=OptActStateSeq( deltas,T,k,u,eta,lambda,gamma,alpha,beta); 
    finalprop(TIdx) = nansum(WithoptActSeqMatrix(:,TIdx));
    U_task_with(TIdx)=alpha*(k*nansum(WithoptActSeqMatrix(:,TIdx)))^beta;
    cost_with(TIdx)=nansum(eta*WithoptActSeqMatrix(:,TIdx).^lambda);
    Utotal_with(TIdx) = U_task_with(TIdx)-cost_with(TIdx);
end

save('T.mat')
% %%%%%%% plot effort as a function of time in different discount rate conditions
% figure(1)
% bluecolor = [158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
% greycolor =[189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
% greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 
% 
% for TIdx=1:length(TVec)
%     T = TVec(TIdx);
%     plot((1:T)',WithoptActSeqMatrix(1:T,TIdx),'o-','Color',greencolor(TIdx,:));
%     hold on
%     legendInfo{TIdx} = ['T=',num2str(T)];
% end
% 
% legend(legendInfo)
% xlabel('time t')
% ylabel('effort a')
% xlim([1,TVec(end)])
% set(gca, 'XTick', 1:TVec(end));
% set(gca, 'FontSize',28)
% 
% %%%%%% plot the U_net, R(s_T+1), total cost
% figure(2)
% plot(TVec,U_task_with,'o-','Color',bluecolor(3,:));
% hold on 
% plot(TVec,cost_with,'o-','Color',redcolor(3,:));
% hold on
% plot(TVec,Utotal_with,'o-','Color',greycolor(3,:));
% 
% legend('reward','total cost','net utility')
% xlabel('total time T')
% ylabel('utility')
% legend boxoff
% set(gca, 'FontSize',28)
% 
% %%%%%%% days of procrastination over alpha  
% figure(3)
% numberzero = nan(1,length(TVec));
% for TIdx = 1:length(TVec)
%     numberzero(TIdx)= nnz(~WithoptActSeqMatrix(1:TVec(TIdx),TIdx))/TVec(TIdx);
% end
% plot(TVec,numberzero,'ko-')
% xlabel('total time T')
% ylabel('normalized days of procrastination')
% ylim([0 1])
% set(gca,'FontSize',28);