%%%%%%%%% switch from not working at all to completing the task in several
%%%%%%%%% senarios
close all;
clear
%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);
figure

%%%%%%%% alpha
deltas = 0.01; k=1;u=0;
gamma=0.6; lambda=3; T=10;
c_1 = 10;
alphaVec = 0:0.1:1;
betaVec = 0.5:0.5:3.5; 
plotbetaIndex = [1,4,7];

WithoptActSeqMatrix=nan(T,length(alphaVec),length(betaVec));
meanActSeqMatrix = nan(length(betaVec),length(alphaVec));
for betaIdx=1:length(betaVec)
    for alphaIdx=1:length(alphaVec)
            
    [WithoptActSeqMatrix(:,alphaIdx,betaIdx), ]=OptActStateSeq( deltas,T,k,u,c_1,lambda,gamma,alphaVec(alphaIdx),betaVec(betaIdx));  
    meanActSeqMatrix(betaIdx,alphaIdx) = nanmean(WithoptActSeqMatrix(:,alphaIdx,betaIdx));
    end
end
greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

figure
plot(alphaVec,k*T*meanActSeqMatrix(plotbetaIndex(1),:),'-o','Color',greencolor(2,:));
hold on
plot(alphaVec,k*T*meanActSeqMatrix(plotbetaIndex(2),:),'o-','Color',greencolor(4,:));
hold on
plot(alphaVec,k*T*meanActSeqMatrix(plotbetaIndex(3),:),'o-','Color',greencolor(6,:));
hold off

legend('\beta=0.5','\beta=2','\beta=3.5')
legend boxoff
xlabel('maximum task reward \alpha')
ylabel('final proportion completed')
xlim([alphaVec(1),alphaVec(end)])
set(gca, 'XTick', alphaVec(1:2:end));
axis square
box off;
% plot(alphaVec,finalprop,'ko-');
% xlabel('maximum task reward \it\alpha','Interpreter','tex')
% ylabel('final proportion completed')
% xlim([alphaVec(1),alphaVec(end)])
% ylim([0,1])
% xticks(alphaVec(1):0.2:alphaVec(end))
% yticks(0:0.2:1)

% %%%%%%% T
% deltas = 0.01; k=1;u=0;
% gamma=0.6; lambda=3;
% beta=4; 
% c_1 = 10;alpha=1;
% TVec = 1:1:10;
% 
% WithoptActSeqMatrix=nan(T,length(TVec));
% finalprop = nan(1,length(TVec));
% 
% for Idx=1:length(TVec)
%     [WithoptActSeqMatrix(1:Idx,Idx), ]=OptActStateSeq( deltas,Idx,k,u,c_1,lambda,gamma,alpha,beta);  
%     finalprop(Idx) = nansum(WithoptActSeqMatrix(1:Idx,Idx));
% end
% 
% subplot(2,3,2)
% plot(TVec,finalprop,'ko-');
% xlabel('total time \itT','Interpreter','tex')
% ylabel('final proportion completed')
% xlim([1,TVec(end)])
% ylim([0,1])
% xticks(TVec)
% yticks(0:0.2:1)
% 
% %%%%%%% lambda
% deltas = 0.01; k=1;u=0;
% gamma=0.6; 
% beta=4; T=10;alpha = 1;
% c_1 = 10;
% lambdaVec = 1.2:0.4:5.2;
% 
% WithoptActSeqMatrix=nan(T,length(lambdaVec));
% finalprop = nan(1,length(lambdaVec));
% 
% for Idx=1:length(lambdaVec)
%     [WithoptActSeqMatrix(1:T,Idx), ]=OptActStateSeq( deltas,T,k,u,c_1,lambdaVec(Idx),gamma,alpha,beta);  
%     finalprop(Idx) = nansum(WithoptActSeqMatrix(1:T,Idx));
% end
% 
% subplot(2,3,3)
% plot(lambdaVec,finalprop,'ko-');
% xlabel('exponent of cost function \it\lambda','Interpreter','tex')
% ylabel('final proportion completed')
% xlim([lambdaVec(1),lambdaVec(end)])
% ylim([0,1])
% xticks(lambdaVec(1):0.8:lambdaVec(end))
% yticks(0:0.2:1)
% 
% %%%%%%% Cmax
% deltas = 0.01; k=1;u=0;
% gamma=0.6; lambda = 2;
% beta=4; T=10;alpha = 1;
% c_1Vec = 1:0.5:6;
% 
% WithoptActSeqMatrix=nan(T,length(c_1Vec));
% finalprop = nan(1,length(c_1Vec));
% 
% for Idx=1:length(c_1Vec)
%     [WithoptActSeqMatrix(1:T,Idx), ]=OptActStateSeq( deltas,T,k,u,c_1Vec(Idx),lambda,gamma,alpha,beta);  
%     finalprop(Idx) = nansum(WithoptActSeqMatrix(1:T,Idx));
% end
% 
% subplot(2,3,4)
% plot(c_1Vec,finalprop,'ko-');
% xlabel('maximum cost \itC_{max}','Interpreter','tex')
% ylabel('final proportion completed')
% xlim([c_1Vec(1),c_1Vec(end)])
% ylim([0,1])
% xticks(c_1Vec(1):1:c_1Vec(end))
% yticks(0:0.2:1)
% 
% %%%%%%%% J
% deltas = 0.01; k=1;
% gamma=0.7; lambda = 3;
% beta=4; T=10;alpha = 1; c_1=7;
% Jvec = 0:0.01:0.08;
% 
% WithoptActSeqMatrix=nan(T,length(Jvec));
% finalprop = nan(1,length(Jvec));
% 
% for Idx=1:length(Jvec)
%     [WithoptActSeqMatrix(1:T,Idx), ]=OptActStateSeq( deltas,T,k,Jvec(Idx),c_1,lambda,gamma,alpha,beta);  
%     finalprop(Idx) = nansum(WithoptActSeqMatrix(1:T,Idx));
% end
% 
% subplot(2,3,5)
% plot(Jvec,finalprop,'ko-');
% xlabel('utility of alternative activities \itJ','Interpreter','tex')
% ylabel('final proportion completed')
% xlim([Jvec(1),Jvec(end)])
% ylim([0,1])
% xticks(Jvec(1):0.02:Jvec(end))
% yticks(0:0.2:1)
% 
% %%%%%% switch from delayed reward to instantaneous reward
% deltas = 0.01; k=1;u=0;
% gamma=0.4; lambda = 3;
% beta=4; T=5;alpha = 1; c_1=7;
% 
% finalprop = nan(1,2);
% 
% [tempdelayed, ]=OptActStateSeq( deltas,T,k,u,c_1,lambda,gamma,alpha,beta);  
% finalprop(1) = nansum(tempdelayed);
% [tempinstan, ]=OptActStateSeqNoDelay( deltas,T,k,u,c_1,lambda,gamma,alpha,beta);  
% finalprop(2) = nansum(tempinstan);
% 
% subplot(2,3,6)
% plot([1,3],finalprop,'ko-')
% ylabel('final proportion completed','Interpreter','tex')
% xlim([0,4])
% xticks([1,3])
% xticklabels({'delayed reward','instantaneous reward'})
% ylim([0,1])
% yticks(0:0.2:1)

