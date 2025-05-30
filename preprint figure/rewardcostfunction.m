
close all;

%% Prepare figures

% Default settings
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

greycolor = [217,217,217;150,150,150;37,37,37; 0,0,0]/255;
%% figure1a reward function 
figure
subplot(1,2,1)
effort = 0:0.001:1;rewardexp=[0.6,3,10];alpha=[0.9,0.8,0.7];
reward=nan(length(rewardexp),length(effort));
reward(1,:)=alpha(1)*effort.^rewardexp(1);
reward(2,:)=alpha(2)*effort.^rewardexp(2);
reward(3,:)=alpha(3)*effort.^rewardexp(3);
%reward(4,:)=alpha(4)*effort.^rewardexp(4);
for iplot=1:length(rewardexp)
plot(effort, reward(iplot,:), '-','Color',greycolor(iplot,:))
hold on
end
legend(['\beta=',num2str(rewardexp(1))],['\beta=',num2str(rewardexp(2))],['\beta=',num2str(rewardexp(3))])%,['\beta=',num2str(rewardexp(4))])
axis([0,1,0,1])
axis square
xlabel('proportion completed \its','Interpreter','tex'); ylabel('task utility \itU','Interpreter','tex')
set(gca,'xtick',[0,1],'ytick',[0,1])

%% figure1b cost function
subplot(1,2,2)
%figure
lambda=0.3;
cost = 0.9*effort.^lambda;
plot(effort, cost, '-','Color',[216,179,101]/255)
hold on
plot([0,1],[0,0.9],'k--')
hold on
%lambda=3;
lambda=4;
cost = 0.9*effort.^lambda;
plot(effort,cost,'-','Color',[90,180,172]/255)
axis([0,1,0,1])
axis square
xlabel('effort \ita','Interpreter','tex'); ylabel('cost \itC','Interpreter','tex')
text(0.1,0.8,'concave','Color',[216,179,101]/255,'FontSize',20)
text(0.6,0.1,'convex','Color',[90,180,172]/255,'FontSize',20)
set(gca,'xtick',[0,1],'ytick',0)
%set(gca,'xtick',[0,0.5,1],'ytick',0)

