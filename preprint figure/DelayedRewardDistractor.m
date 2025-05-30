deltas = 0.01; k=1;uvec=0:0.002:0.024;
eta=7; lambda=3;
alpha=1; beta=2;
gamma = 0.6; T = 10;

WithoptActSeqMatrix=nan(T,length(uvec));
meanActSeqMatrix = nan(1,length(uvec));
NormDaysofPro = nan(1,length(uvec));
U_task_with=nan(1,length(uvec));
cost_with=nan(1,length(uvec));
Utotal_with=nan(1,length(uvec));
finalprop = nan(1,length(uvec));

for uIdx = 1:length(uvec)
    WithoptActSeqMatrix(:,uIdx)=OptActStateSeq( deltas,T,k,uvec(uIdx),eta,lambda,gamma,alpha,beta);
    U_task_with(1,uIdx)=alpha*(nansum(WithoptActSeqMatrix(:,uIdx)))^beta;
    cost_with(1,uIdx)=nansum(eta*WithoptActSeqMatrix(:,uIdx).^lambda);
    Utotal_with(1,uIdx) = U_task_with(:,uIdx)-cost_with(:,uIdx);
    finalprop(1,uIdx) = nansum(WithoptActSeqMatrix(:,uIdx));
end


save('J.mat')
%save('distrtimecoursecurve_daysofpro_utility.mat')
% 
% figure
% bluecolor = [158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
% greycolor =[189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0]/255;
% greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
% redcolor = [252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]/255; 
% 
% plotuIdx = [1,2,3,4,5];
% for uIdx = 1:length(plotuIdx)
%     plot(1:T,WithoptActSeqMatrix(:,plotuIdx(uIdx)),'o-','Color',greencolor(uIdx,:));
%     hold on
%     legendInfo{uIdx} = ['r_{dis}=',num2str(uvec(plotuIdx(uIdx)))];
% end
% xlim([1,T])
% legend(legendInfo)
% xlabel('time t')
% ylabel('effort a')
% set(gca, 'XTick', 1:T,'FontSize',28);
% 
% figure
% numberzero = nan(1,length(uvec));
% for inumber = 1:length(uvec)
% numberzero(inumber)=nnz(~WithoptActSeqMatrix(:,inumber));
% end
% plot(uvec,numberzero,'ko-')
% xlabel('distractor reward r_{dis}')
% ylabel('days of procrastination')
% xlim([min(uvec),max(uvec)])
% set(gca,'FontSize',28);
% 
% 
% figure
% plot(uvec,U_task_with,'o-','Color',bluecolor(3,:));
%     hold on 
%     plot(uvec,cost_with,'o-','Color',redcolor(3,:));
%     hold on
%     plot(uvec,Utotal_with,'o-','Color',greycolor(3,:));
% 
% legend('reward','total cost','net utility')
% xlabel('distractor reward r_{dis}')
% ylabel('utility')
% xlim([min(uvec),max(uvec)])
% set(gca, 'FontSize',28,'XTick',uvec)
% 
