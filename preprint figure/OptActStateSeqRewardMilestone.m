function [ OptActSequence,OptStateSequence,Value ] = OptActStateSeqRewardMilestone( Deltas,T,k,r_0,p,q,gamma,m,n,N_milestone )
% this version is revised based on
% OptActStateSeqRewardProgress_Incorrect_NeedtobeRevised.m 
% the change I make is separate T from before T, because in the wrong
% version, if an agent did not complete the task, their additional progress
% beyond the latest milestone will not be rewarded. But in this version, I
% corrected. Additionally, I added N_milestone, previously N_milestone=T (a special case) in OptActStateSeqRewardProgress_Incorrect_NeedtobeRevised.m 

% this is the version without delayed feedback U_t = alpha*s_t^beta  
% deltas = k*a
% Deltas      =  0.01; % should be 0.01
StateVector =  (0:Deltas:1)';
% T           =  5;
StateMatrix =  nan(length(StateVector),T+1);
StateMatrix(1,1)=0;
StateMatrix(:,2:end)=repmat(StateVector,1,T);

% action deltas=k*a
% k = 1; 
ActVector = StateVector/k; % because state unit is fixed here, so k<=1,
% if k>1, by having an action ActVector, the state will end up somewhere
% within the unit. 

% reward values
% r_0 = 0;
% R = s, could be a separate function, but for simplicity here, R=s. 
% p=2; q=6;
cost = p*ActVector.^q; 
% cost =  p*ActVector; 
% cost(q:end)=2;
% myopia parameter
% m =1; n=1;
% gamma = 0.2;

% Let's say reward every 1/T = 1/5 
%reward_interval_list = [0.25,0.5,0.75,1];
reward_interval_list = cumsum(ones(1,N_milestone)*1/N_milestone);
%reward_interval_list = [0.1, 0.3, 0.9, 1];

Value       =nan(length(StateVector),T+1);
Value(:,T+1)=0;
OptActIn    =nan(length(StateVector),T); % choose the action with index in ActVector.
%% First calculate Q(s,a) at T
% on the last day, if the agent did not complete the task, they will be rewarded
% according to their additional progress from the last rewarding time point. (if the agent complete the task, then have fun instead)

% calculate the current state<1 (not complete the task)
for stIn = 1:size(StateMatrix,1)-1 % if the current state not yet complete the task, if complete the task stIn= size(StateMatrix,1)
    NextStateVector = StateMatrix(stIn,T)+k*ActVector;
    AvailActIndex = find(NextStateVector<=1);
    Q = nan(length(AvailActIndex),1);
    
    % if have fun
    % additional reward = reward for final s - reward for last rewarded s
    % now let's calculate reward for last rewarded s
    listminusstate = reward_interval_list - StateMatrix(stIn,T);
    closestindexOnLeftOfState = sum(listminusstate<=0);
    if closestindexOnLeftOfState == 0 % meaning StateMatrix(stIn,T)<reward_interval_list(1)
        U_Ts = 0;
    else
        U_Ts = m*reward_interval_list(closestindexOnLeftOfState)^n;
    end
    
    % now let's calculate reward for final s
    U_Tplus1 = m*StateMatrix(stIn,T+1)^n;
    if U_Tplus1>1
        U_Tplus1=1; % U is constrained between 0 and 1
    end
    
    % additional reward = reward for final s - reward for last rewarded s
    Q(1)= r_0+U_Tplus1-U_Ts+gamma*Value(stIn,T+1);
    
    % if work 
    for QIn = 2: size(Q,1)
        U_tplus1sprime = m*StateMatrix(stIn+QIn-1,T+1)^n;
        if U_tplus1sprime>1
            U_tplus1sprime=1; % U is constrained between 0 and 1
        end
        % additional reward = reward for final s - reward for last rewarded s
        Q(QIn)= -cost(QIn)+U_tplus1sprime-U_Ts+gamma*Value(stIn+QIn-1,T+1);
    end
    [Value(stIn,T),OptActIn(stIn,T)]=max(Q);    
end
% if current state is task completed:   
OptActIn(size(StateMatrix,1),T) = 1; % choose for fun (only choice when you already completed the task at T)
Value(size(StateMatrix,1),T) = r_0+gamma*Value(size(StateMatrix,1),T+1);

for t = T-1:-1:1
    if t==1
        stInend = 1;
    else 
        stInend = size(StateMatrix,1);
    end
    for stIn = 1:stInend
        NextStateVector = StateMatrix(stIn,t)+k*ActVector;
        AvailActIndex = find(NextStateVector<=1);
        Q = nan(length(AvailActIndex),1);
        % for fun ActVector(1)=0;
        %U_tplus1s = m*StateMatrix(stIn,t+1)^n;
        
        %U_ts = m*StateMatrix(stIn,t)^n;
        
        % mathematically U_tplus1s always = U_ts
        Q(1)= r_0+gamma*Value(stIn,t+1); % Q(1)= r_0+U_tplus1s-U_ts+gamma*Value(stIn,t+1); 
        
        % for work
        listminusstate = reward_interval_list - StateMatrix(stIn,t);
        closestindexOnLeftOfState = sum(listminusstate<=0); % listminusstate<=0 give us logic [1,1,0,0,0] for example, sum of listminusstate<=0 gives us index of the closeset element in reward_interval_list to the left side of StateMatrix(stIn,t)
        if closestindexOnLeftOfState == 0 % meaning StateMatrix(stIn,t)<reward_interval_list(1)
            U_ts = 0;
        else
            U_ts = m*reward_interval_list(closestindexOnLeftOfState)^n;
        end
% the following code is used before I know how to code in a more efficient way as the above        
%         if StateMatrix(stIn,t)<reward_interval_list(1)
%             U_ts = 0;
%         end
%         if StateMatrix(stIn,t)>=reward_interval_list(1) && StateMatrix(stIn,t)<reward_interval_list(2)
%             U_ts = m*reward_interval_list(1)^n;
%         end
%         if StateMatrix(stIn,t)>=reward_interval_list(2) && StateMatrix(stIn,t)<reward_interval_list(3)
%             U_ts = m*reward_interval_list(2)^n;
%         end        
%         if StateMatrix(stIn,t)>=reward_interval_list(3) && StateMatrix(stIn,t)<reward_interval_list(4)
%             U_ts = m*reward_interval_list(3)^n;
%         end
%         if StateMatrix(stIn,t)>=reward_interval_list(4) && StateMatrix(stIn,t)<reward_interval_list(5)
%             U_ts = m*reward_interval_list(4)^n;
%         end      
%         if StateMatrix(stIn,t)==reward_interval_list(end)
%             U_ts = m*reward_interval_list(end)^n;
%         end
        
        for QIn = 2: size(Q,1)
            listminusstate2 = reward_interval_list - StateMatrix(stIn+QIn-1,t+1);
            closestindexOnLeftOfState2 = sum(listminusstate2<=0);
            if closestindexOnLeftOfState2 == 0
                U_tplus1sprime = 0;
            else
                U_tplus1sprime = m*reward_interval_list(closestindexOnLeftOfState2)^n;
            end
% the following code is used before I know how to code in a more efficient way as the above       
%             if StateMatrix(stIn+QIn-1,t+1)<reward_interval_list(1)
%                 U_tplus1sprime = 0;
%             end            
%             if StateMatrix(stIn+QIn-1,t+1)>=reward_interval_list(1) && StateMatrix(stIn+QIn-1,t+1)<reward_interval_list(2)
%                 U_tplus1sprime = m*reward_interval_list(1)^n;
%             end
%             if StateMatrix(stIn+QIn-1,t+1)>=reward_interval_list(2) && StateMatrix(stIn+QIn-1,t+1)<reward_interval_list(3)
%                 U_tplus1sprime = m*reward_interval_list(2)^n;
%             end        
%             if StateMatrix(stIn+QIn-1,t+1)>=reward_interval_list(3) && StateMatrix(stIn+QIn-1,t+1)<reward_interval_list(4)
%                 U_tplus1sprime = m*reward_interval_list(3)^n;
%             end
%             if StateMatrix(stIn+QIn-1,t+1)>=reward_interval_list(4) && StateMatrix(stIn+QIn-1,t+1)<reward_interval_list(5)
%                 U_tplus1sprime = m*reward_interval_list(4)^n;
%             end 
%             if StateMatrix(stIn+QIn-1,t+1) == reward_interval_list(end)
%                 U_tplus1sprime = m*reward_interval_list(end)^n;
%             end
            %U_tplus1sprime = m*StateMatrix(stIn+QIn-1,t+1)^n;
%             if U_tplus1sprime>1
%                 U_tplus1sprime=1; % U is constrained between 0 and 1
%             end
            
            Q(QIn)= U_tplus1sprime-U_ts-cost(QIn)+gamma*Value(stIn+QIn-1,t+1);
        end
        [Value(stIn,t),OptActIn(stIn,t)]=max(Q);
        %Value(stIn,t)=max(Q);
        %tempOptActIn = find(abs(Q - Value(stIn,t))<eps);
        %tempSeq = randperm(length(tempOptActIn));
        %OptActIn(stIn,t) = tempSeq(1);
    end
end

%% find the optimal action sequence and its corresponding state sequence
% with state action trainsition rule 
OptStateSequence = nan(1,T+1);
OptStateSequence(1)=0;
OptStateSequence(2) = ActVector(OptActIn(1,1))*k; 

OptActSeqenceIndex=nan(1,T);
OptActSeqenceIndex(1)=OptActIn(1,1);
for t = 2: T
    indexTemp=find(abs(StateMatrix(:,t)-OptStateSequence(t))<0.000001);
    
    OptActSeqenceIndex(t)=OptActIn(indexTemp,t);
    delta_s= ActVector(OptActSeqenceIndex(t))*k;
    OptStateSequence(t+1)=OptStateSequence(t)+delta_s;
end
OptActSequence = ActVector(OptActSeqenceIndex); 


end

