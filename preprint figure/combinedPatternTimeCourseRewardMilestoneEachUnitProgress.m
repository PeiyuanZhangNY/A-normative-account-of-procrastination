%% check the total number of patterns of work and figure for patterns of time course of progress in reward-each-milestone schedule and reward each unit of completion 

load('ImmediateRewardIntervention_2.mat')

%% check the total number of patterns of work

[labels, counts, names,patternMatrices] = classifyTimeCoursePatterns(WithoptActSeqMatrix3);
% Pattern classification of 10000 time courses (T = 10 days):
% Pattern           Count  Percent
% Steady             1527    15.3%
% RampUp              503     5.0%
% RampDown           5775    57.8%
% UpThenDown         1055    10.6%
% DownThenUp          458     4.6%
% Fluctuating         682     6.8%
% Distinct pattern kinds present: 6

[labels, counts, names,patternMatrices] = classifyTimeCoursePatterns(WithoptActSeqMatrix4);
% Pattern classification of 10000 time courses (T = 10 days):
% Pattern           Count  Percent
% Steady             1485    14.8%
% RampUp               14     0.1%
% RampDown           7432    74.3%
% UpThenDown         1069    10.7%
% DownThenUp            0     0.0%
% Fluctuating           0     0.0%
% Distinct pattern kinds present: 4

%% figure
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);

greencolor = [161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;

figure
subplot(2,5,1)
T=10;
% working steadily
i=8;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
plot(1:T,plotVector,'o-','Color',greencolor(4,:));
%xlabel('time \itt','Interpreter','tex')
%ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.8])
xticks(1:T)
yticks(0:0.2:0.8)
set(gca,'TickDir','out');
%axis square

subplot(2,5,6)
% ramping up 
i=14;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
plot(1:T,plotVector,'o-','Color',greencolor(4,:)); 
%xlabel('time \itt','Interpreter','tex')
%ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.8])
xticks(1:T)
yticks(0:0.2:0.8)
set(gca,'TickDir','out');
%axis square

subplot(2,5,2)
% going down
i = 7;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
plot(1:T,plotVector,'o-','Color',greencolor(4,:));  
%xlabel('time \itt','Interpreter','tex')
%ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.8])
xticks(1:T)
yticks(0:0.2:0.8)
set(gca,'TickDir','out');
%axis square

subplot(2,5,7)
% first going down and then ramping up (U=shape)
i = 49;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
plot(1:T,plotVector,'o-','Color',greencolor(4,:));
%xlabel('time \itt','Interpreter','tex')
%ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.8])
xticks(1:T)
yticks(0:0.2:0.8)
set(gca,'TickDir','out');
%axis square

subplot(2,5,3)
% first ramping up and then going down (inverted U=shape)
i = 74;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
plot(1:T,plotVector,'o-','Color',greencolor(4,:));
%xlabel('time \itt','Interpreter','tex')
%ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.8])
xticks(1:T)
yticks(0:0.2:0.8)
set(gca,'TickDir','out');
%axis square

subplot(2,5,8)
% fluctuate
i = 87;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,3)) = WithoptActSeqMatrix3(i,1:finisheday(i,3));
plot(1:T,plotVector,'o-','Color',greencolor(4,:)); 
xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.3])
xticks(1:T)
yticks(0:0.1:0.3)
set(gca,'TickDir','out');
%axis square


%% figure for patterns of time course of progress in reward-each-unit-progress schedule 

subplot(2,5,4)
T=10;
% working steadily
i=63;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,4)) = WithoptActSeqMatrix4(i,1:finisheday(i,4));
plot(1:T,plotVector,'o-','Color',greencolor(2,:));
%xlabel('time \itt','Interpreter','tex')
ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.6])
xticks(1:T)
yticks(0:0.2:0.6)
set(gca,'TickDir','out');
%axis square

subplot(2,5,9)
% ramping up 
i=3;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,4)) = WithoptActSeqMatrix4(i,1:finisheday(i,4));
plot(1:T,plotVector,'o-','Color',greencolor(2,:));  
%xlabel('time \itt','Interpreter','tex')
%ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.6])
xticks(1:T)
yticks(0:0.2:0.6)
set(gca,'TickDir','out');
%axis square

subplot(2,5,5)
% going down
i = 13;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,4)) = WithoptActSeqMatrix4(i,1:finisheday(i,4));
plot(1:T,plotVector,'o-','Color',greencolor(2,:));  
%xlabel('time \itt','Interpreter','tex')
%ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.6])
xticks(1:T)
yticks(0:0.2:0.6)
set(gca,'TickDir','out');
%axis square

subplot(2,5,10)
% first ramping up and then going down (inverted U=shape)
i = 5;
plotVector = NaN(1, T); 
plotVector(1:finisheday(i,4)) = WithoptActSeqMatrix4(i,1:finisheday(i,4));
plot(1:T,plotVector,'o-','Color',greencolor(2,:));  
xlabel('time \itt','Interpreter','tex')
%ylabel('progress \Delta\its','Interpreter','tex')
xlim([1,T])
ylim([0,0.6])
xticks(1:T)
yticks(0:0.2:0.6)
set(gca,'TickDir','out');
%axis square

