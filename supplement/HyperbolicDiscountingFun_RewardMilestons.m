function OptActSequence = HyperbolicDiscountingFun_RewardMilestons(T,k,r_0,p,q,m,n,k_DR,N_milestone)

DeltaGammaVector = 0.2;
GammaVector = 0:DeltaGammaVector:0.99; % the maximum has to be 0.99, if it is 1, 
% then even if k_DR=100, then there will be no effort=0's.  because
% 1^100=1, but 0.99^100=0.366. They are much away from each other. And also
% for math correction, should be 0.99 if I use DeltaGammaVector. (left) Riemann sum
% In paper William Fedus et al,2019 
TransGammaVector = GammaVector.^k_DR; % transferred gamma vector 
Deltas = 0.01;StateVector =  (0:Deltas:1)';ActVector = StateVector/k;
StateMatrix =  nan(length(StateVector),T+1);
StateMatrix(1,1)=0;
StateMatrix(:,2:end)=repmat(StateVector,1,T);

%% calculate Q_sa, here, very similar with function OptActStateSeqRewardMilestone with a few modifications marked
reward_interval_list = cumsum(ones(1,N_milestone)*1/N_milestone);
cost = p*ActVector.^q; 
gammaVec = TransGammaVector;
Q_sa = nan(length(gammaVec),T,length(StateVector),length(ActVector));

for gammaIdx = 1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    
    %% First calculate Q(s,a) at T 
    Value       =nan(length(StateVector),T+1);
    Value(:,T+1)=0;    
    
    % calculate the current state<1 (not the final state)
    for stIn = 1:size(StateMatrix,1)-1 % The last state is 1, then it should have been rewarded at T-1, so should be separately analyzed
    NextStateVector = StateMatrix(stIn,T)+k*ActVector;
    AvailActIndex = find(NextStateVector<=1);
 
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
    Q_sa(gammaIdx,T,stIn,1) = r_0+U_Tplus1-U_Ts+gamma*Value(stIn,T+1);

    % if work 
    for ActIn = 2: length(AvailActIndex)
        U_tplus1sprime = m*StateMatrix(stIn+ActIn-1,T+1)^n;
        if U_tplus1sprime>1
            U_tplus1sprime=1; % U is constrained between 0 and 1
        end
        % additional reward = reward for final s - reward for last rewarded s
        Q_sa(gammaIdx,T,stIn,ActIn)= -cost(ActIn)+U_tplus1sprime-U_Ts+gamma*Value(stIn+ActIn-1,T+1);
    end
    Value(stIn,T)=nanmax(Q_sa(gammaIdx,T,stIn,:));
    end
    
    % if current state is task completed:   
    Value(size(StateMatrix,1),T) = r_0+gamma*Value(size(StateMatrix,1),T+1);

    %% Second calculate Q(s,a) at t<T
    for t = T-1:-1:1
    if t==1
        stInend = 1;
    else 
        stInend = size(StateMatrix,1);
    end   
    
    for stIn = 1:stInend
        NextStateVector = StateMatrix(stIn,t)+k*ActVector;
        AvailActIndex = find(NextStateVector<=1);
        % for fun ActVector(1)=0;
        %U_tplus1s = m*StateMatrix(stIn,t+1)^n;
        
        %U_ts = m*StateMatrix(stIn,t)^n;
        
        % mathematically U_tplus1s always = U_ts
        Q_sa(gammaIdx,t,stIn,1)= r_0+gamma*Value(stIn,t+1); % Q(1)= r_0+U_tplus1s-U_ts+gamma*Value(stIn,t+1); 
        
        % for work
        listminusstate = reward_interval_list - StateMatrix(stIn,t);
        closestindexOnLeftOfState = sum(listminusstate<=0); % listminusstate<=0 give us logic [1,1,0,0,0] for example, sum of listminusstate<=0 gives us index of the closeset element in reward_interval_list to the left side of StateMatrix(stIn,t)
        if closestindexOnLeftOfState == 0 % meaning StateMatrix(stIn,t)<reward_interval_list(1)
            U_ts = 0;
        else
            U_ts = m*reward_interval_list(closestindexOnLeftOfState)^n;
        end
        
        for ActIn = 2: length(AvailActIndex)
            listminusstate2 = reward_interval_list - StateMatrix(stIn+ActIn-1,t+1);
            closestindexOnLeftOfState2 = sum(listminusstate2<=0);
            if closestindexOnLeftOfState2 == 0
                U_tplus1sprime = 0;
            else
                U_tplus1sprime = m*reward_interval_list(closestindexOnLeftOfState2)^n;
            end
            Q_sa(gammaIdx,t,stIn,ActIn)=U_tplus1sprime-U_ts-cost(ActIn)+gamma*Value(stIn+ActIn-1,t+1);
        end
        Value(stIn,t)=max(Q_sa(gammaIdx,t,stIn,:));
    end
    end
end

%% calculate Q_saK
Q_saK = nan(T,length(StateVector),length(ActVector));
OptActIn = nan(length(StateVector),T);

for t = 1:T
    for stIn = 1:size(StateMatrix,1)
        NextStateVector = StateMatrix(stIn,T)+round(k*ActVector/Deltas)*Deltas;
        AvailActIndex = find(NextStateVector<=1);
        
        for ActIn = 1: length(AvailActIndex)
            Q_saK(t,stIn,ActIn) = nansum(Q_sa(:,t,stIn,ActIn))*DeltaGammaVector;
        end
        [~,OptActIn(stIn,t)] = nanmax(Q_saK(t,stIn,:));
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