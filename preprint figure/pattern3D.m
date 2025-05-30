%%% the original file is in TheoryPaper1_2Figures folder pattern3D.m, and
%%% pattern3D.pdf

%3D of all the possible behavior pattern including 
% 1) work last-minute, 2) not working at all, 3) a delay in the beginning and rushing in the end, 4) increasing trend
% parameters lambda, gamma, c1
%%%%% in delayed reward condition
% when lambda<1, work last-minute,alpha>c1
% when lambda>1, alpha<c1, not working at all
% when lambda>1, gamma small, a delay...
% when lambda>1, gamma large, increasing trend
close all; 
clear;

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);
%% in delayed reward condition
start=0.1;deltapara=0.3; last =2.6;
lambdaVec = start:deltapara:last; % 0.1:0.4:5;
c1Vec = lambdaVec;
%gammaVec = [0.1:0.3:0.7,0.99];
gamma = 0.7;
alpha = 1; beta=1.9; J=0; k=1; deltas=0.01;T = 6;
patternDelay = nan(length(lambdaVec),length(c1Vec));
for idxlambda = 1:length(lambdaVec)
    for idxc1 = 1:length(c1Vec)
        %for idxgamma = 1:length(gammaVec)
            [temp,] = OptActStateSeq(deltas,T,k,J,c1Vec(idxc1),lambdaVec(idxlambda),gamma,alpha,beta); 
            if sum(temp)==0 
                patternDelay(idxlambda,idxc1) = 0; % not working at all
                
            elseif find(temp,1)==T
                patternDelay(idxlambda,idxc1) = 0.5; % work last-minute

            elseif find(temp,1)<T && find(temp,1)>1  
                patternDelay(idxlambda,idxc1) = 1; % a delay ...
    
            else
                patternDelay(idxlambda,idxc1) = 1.5; % no delay 
            end
            
            
        %end
    end
end

%graph1 = squeeze(patternDelay(:,:,1));
%graph2 = squeeze(patternDelay(:,:,2));
graph2 = patternDelay;
%graph3 = squeeze(patternDelay(:,:,3));

%% imagesec for different gamma
figure 
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);
% %map = [215,25,28;171,217,233;44,123,182;82,168,38;255,255,191]/255;
% reddarklightblue = [215,25,28;44,123,182;171,217,233]/255; % red:215,25,28;(not working at all) light blue:171,217,233; (work last minute) dark blue: 44,123,182 (a delay)  
% ax(1) = subplot(1,3,1);
% yaxislambda = lambdaVec(1)-deltapara/2:deltapara/2:lambdaVec(end)-deltapara/2;
% xaxisc1 = c1Vec(1)-deltapara/2:deltapara/2:c1Vec(end)-deltapara/2;
% imagesc(xaxisc1, yaxislambda, graph1);
% colormap(ax(1),reddarklightblue)
% set(gca, 'XTick', c1Vec,'YTick',lambdaVec);
% grid on
% ax = gca;
% ax.GridColor = [0 0 0];
% ax.GridAlpha = 1;
% ax.LineWidth = 2;
% axis square;
% ylabel('\it\lambda','Interpreter','tex'); xlabel('\itc_1','Interpreter','tex')
% title('\gamma=0.1','Interpreter','tex')

% ax(2) = subplot(1,3,2);
yaxislambda = lambdaVec(1)-deltapara/2:deltapara/2:lambdaVec(end)-deltapara/2;
xaxisc1 = c1Vec(1)-deltapara/2:deltapara/2:c1Vec(end)-deltapara/2;
imagesc(xaxisc1, yaxislambda, graph2);
reddarklightblueyellow = [215,25,28;253,174,97;171,217,233;44,123,182]/255; % yellow: 255,255,191 (no delay)
colormap(reddarklightblueyellow)
set(gca, 'XTick', c1Vec,'YTick',lambdaVec);
%xlim([0,Tven(end)])
grid on
ax = gca;
ax.GridColor = [0 0 0];
ax.GridAlpha = 1;
ax.LineWidth = 2;
axis square;
ylabel('exponent of cost function \it\lambda','Interpreter','tex'); xlabel('maximum cost \itc_{max}','Interpreter','tex')
%title('\gamma=0.4','Interpreter','tex')
x0=10;
y0=10;
width=550;
height=400;
set(gcf,'position',[x0,y0,width,height])
% ax(3) =subplot(1,3,3);
% yaxislambda = lambdaVec(1)-deltapara/2:deltapara/2:lambdaVec(end)-deltapara/2;
% xaxisc1 = c1Vec(1)-deltapara/2:deltapara/2:c1Vec(end)-deltapara/2;
% imagesc(xaxisc1, yaxislambda, graph3);
% reddarkblueyellow = [215,25,28;44,123,182;255,255,191]/255;
% colormap(ax(3),reddarklightblueyellow)
% set(gca, 'XTick', c1Vec,'YTick',lambdaVec);
% grid on
% ax = gca;
% ax.GridColor = [0 0 0];
% ax.GridAlpha = 1;
% ax.LineWidth = 2;
% axis square;
% ylabel('\it\lambda','Interpreter','tex'); xlabel('\itc_1','Interpreter','tex')
% title('\gamma=0.7','Interpreter','tex')


%% plot the examples for delayed reward condition
figure(2)
subplot(4,1,4) % not working at all
plot(1:T,zeros(1,T),'o-','Color',reddarklightblueyellow(1,:),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',reddarklightblueyellow(1,:))
xlim([1,T])
ylim([0,1])
xlabel('time \it t','Interpreter','tex')
set(gca,'FontSize',16)

subplot(4,1,3) % working at last minute
plot(1:T,[zeros(1,T-1),1],'o-','Color',reddarklightblueyellow(2,:),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',reddarklightblueyellow(2,:))
xlim([1,T])
ylim([0,1])
ylabel('progress \it \Delta s','Interpreter','tex')
set(gca,'FontSize',16)

subplot(4,1,2) % a delay in the beginning
idxlambda = 5;idxc1=1;
[temp,] = OptActStateSeq(deltas,T,k,J,c1Vec(idxc1),lambdaVec(idxlambda),gamma,alpha,beta);
plot(1:T,temp,'o-','Color',reddarklightblueyellow(3,:),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',reddarklightblueyellow(3,:))
xlim([1,T])
ylim([0,1])
set(gca,'FontSize',16)

subplot(4,1,1) % no delay
idxlambda = length(lambdaVec);idxc1=1;
[temp,] = OptActStateSeq(deltas,T,k,J,c1Vec(idxc1),lambdaVec(idxlambda),gamma,alpha,beta);
plot(1:T,temp,'o-','Color',reddarklightblueyellow(4,:),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',reddarklightblueyellow(4,:))
xlim([1,T])
ylim([0,0.4])
set(gca,'FontSize',16)

set(gcf,'position',[x0,y0,width,height])


% %% plot in the instantaneous reward condition
% start=0.1;deltapara=0.3; last =2.6;
% lambdaVec = start:deltapara:last; % 0.1:0.4:5;
% c1Vec = lambdaVec;
% gammaVec = [0.1:0.3:0.7,0.99];
% alpha = 1; beta=1.9; J=0; k=1; deltas=0.01;T = 5;
% patternInstan = nan(length(lambdaVec),length(c1Vec),length(gammaVec));
% for idxlambda = 1:length(lambdaVec)
%     for idxc1 = 1:length(c1Vec)
%         for idxgamma = 1:length(gammaVec)
%             [temp,] = OptActStateSeqNoDelay(deltas,T,k,J,c1Vec(idxc1),lambdaVec(idxlambda),gammaVec(idxgamma),alpha,beta); 
%             if sum(temp)==0 % not working at all
%                 patternInstan(idxlambda,idxc1,idxgamma) = 0;
%                 
%             elseif find(temp==0, 1, 'first')==2
%                 patternInstan(idxlambda,idxc1,idxgamma) = 0.5; % work only on the first day
%    
%             else
%                 patternInstan(idxlambda,idxc1,idxgamma) = 1; % no delay 
%             end
%             
%             
%         end
%     end
% end
% 
% graph4 = squeeze(patternInstan(:,:,1));
% graph5 = squeeze(patternInstan(:,:,2));
% graph6 = squeeze(patternInstan(:,:,3));
% 
% 
% %% imagesec for different gamma
% 
% ax(4)=subplot(2,3,4);
% yaxislambda = lambdaVec(1)-deltapara/2:deltapara/2:lambdaVec(end)-deltapara/2;
% xaxisc1 = c1Vec(1)-deltapara/2:deltapara/2:c1Vec(end)-deltapara/2;
% imagesc(xaxisc1, yaxislambda, graph4);
% redgreenyellow = [215,25,28;82,168,38;255,255,191]/255;
% colormap(ax(4),redgreenyellow)
% set(gca, 'XTick', c1Vec,'YTick',lambdaVec);
% %xlim([0,Tven(end)])
% grid on
% ax = gca;
% ax.GridColor = [0 0 0];
% ax.GridAlpha = 1;
% ax.LineWidth = 2;
% axis square;
% ylabel('\it\lambda','Interpreter','tex'); xlabel('\itc_1','Interpreter','tex')
% title('\gamma=0.1','Interpreter','tex')
% 
% ax(5)=subplot(2,3,5);
% yaxislambda = lambdaVec(1)-deltapara/2:deltapara/2:lambdaVec(end)-deltapara/2;
% xaxisc1 = c1Vec(1)-deltapara/2:deltapara/2:c1Vec(end)-deltapara/2;
% imagesc(xaxisc1, yaxislambda, graph5);
% colormap(ax(5),redgreenyellow)
% set(gca, 'XTick', c1Vec,'YTick',lambdaVec);
% %xlim([0,Tven(end)])
% grid on
% ax = gca;
% ax.GridColor = [0 0 0];
% ax.GridAlpha = 1;
% ax.LineWidth = 2;
% axis square;
% ylabel('\it\lambda','Interpreter','tex'); xlabel('\itc_1','Interpreter','tex')
% title('\gamma=0.4','Interpreter','tex')
% 
% ax(6)=subplot(2,3,6);
% yaxislambda = lambdaVec(1)-deltapara/2:deltapara/2:lambdaVec(end)-deltapara/2;
% xaxisc1 = c1Vec(1)-deltapara/2:deltapara/2:c1Vec(end)-deltapara/2;
% imagesc(xaxisc1, yaxislambda, graph6);
% colormap(ax(6),redgreenyellow)
% set(gca, 'XTick', c1Vec,'YTick',lambdaVec);
% %xlim([0,Tven(end)])
% grid on
% ax = gca;
% ax.GridColor = [0 0 0];
% ax.GridAlpha = 1;
% ax.LineWidth = 2;
% axis square;
% ylabel('\it\lambda','Interpreter','tex'); xlabel('\itc_1','Interpreter','tex')
% title('\gamma=0.7','Interpreter','tex')

