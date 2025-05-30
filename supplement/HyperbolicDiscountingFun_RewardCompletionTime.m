function OptActSequence = HyperbolicDiscountingFun_RewardCompletionTime(T,k,r_0,p,q,m,n,k_DR)

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

%% calculate Q_sa, here, very similar with function OptActStateSeqRewardCompletionTime with a few modifications marked

cost = p*ActVector.^q; 
gammaVec = TransGammaVector;
Q_sa = nan(length(gammaVec),T,length(StateVector),length(ActVector));

for gammaIdx = 1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    
    %% First calculate Q(s,a) at T 
    Value       =nan(length(StateVector),T+1);
    Value(:,T+1)=0;
    % because the reward given at T+1 is different from that given before T+1.

    % calculate the current state<1 (not the final state)
    for stIn = 1:size(StateMatrix,1)-1 % The last state is 1, then it should have been rewarded at T-1, so should be separately analyzed
    NextStateVector = StateMatrix(stIn,T)+k*ActVector;
    AvailActIndex = find(NextStateVector<=1);
    
    % for fun ActVector(1)=0;
    U_Tplus1 = m*StateMatrix(stIn,T+1)^n;
    if U_Tplus1>1
        U_Tplus1=1; % U is constrained between 0 and 1
    end
    Q_sa(gammaIdx,T,stIn,1)= r_0+U_Tplus1+gamma*Value(stIn,T+1);
    
    % for work
    for ActIn = 2: length(AvailActIndex)
        U_tplus1sprime = m*NextStateVector(ActIn)^n; % reward
        if U_tplus1sprime>1
            U_tplus1sprime=1; % U is constrained between 0 and 1
        end
        Q_sa(gammaIdx,T,stIn,ActIn)= -cost(ActIn)+U_tplus1sprime+gamma*Value(stIn+ActIn-1,T+1);

    end 
    Value(stIn,T)=nanmax(Q_sa(gammaIdx,T,stIn,:));
    end
    % if the current state is 1 
    Value(size(StateMatrix,1),T) = r_0+gamma*Value(size(StateMatrix,1),T+1);

    %% Second calculate Q(s,a) at t=T-1:-1:2
    for t = T-1:-1:2
    %if t==1
    %    stInend = 1;
    %else 
    stInend = size(StateMatrix,1);
    %end
    % calculate the current state<1 (not the final state)
    for stIn = 1:stInend-1 
        NextStateVector = StateMatrix(stIn,t)+k*ActVector;
        AvailActIndex = find(NextStateVector<=1);
        % for fun ActVector(1)=0;
        Q_sa(gammaIdx,t,stIn,1)= r_0+gamma*Value(stIn,t+1);         
        % for work        
        for ActIn = 2: length(AvailActIndex)
            U_tplus1sprime = m*StateMatrix(stIn+ActIn-1,t+1)^n;
            if U_tplus1sprime>1
               U_tplus1sprime=1; % U is constrained between 0 and 1
            end
            if StateMatrix(stIn+ActIn-1,t+1) == 1 
                Q_sa(gammaIdx,t,stIn,ActIn)= U_tplus1sprime-cost(ActIn)+gamma*Value(stIn+ActIn-1,t+1); % m is the maximum reward value
            else
                Q_sa(gammaIdx,t,stIn,ActIn)= -cost(ActIn)+gamma*Value(stIn+ActIn-1,t+1);
            end
        end
        Value(stIn,t)=max(Q_sa(gammaIdx,t,stIn,:));
        %Value(stIn,t)=max(Q);
        %tempOptActIn = find(abs(Q - Value(stIn,t))<eps);
        %tempSeq = randperm(length(tempOptActIn));
        %OptActIn(stIn,t) = tempSeq(1);
    end
    
    % if the current state is 1 

    Value(size(StateMatrix,1),t) = r_0+gamma*Value(size(StateMatrix,1),t+1);    
    end
    
    %% last calculate Q(s,a) at t=1
    NextStateVector = StateMatrix(1,1)+k*ActVector;
    AvailActIndex = find(NextStateVector<=1);
    % for fun ActVector(1)=0;
    Q_sa(gammaIdx,1,1,1)= r_0+gamma*Value(1,2); 

    % for work        
    for ActIn = 2: length(AvailActIndex)
        U_tplus1sprime = m*StateMatrix(ActIn,2)^n;
        if U_tplus1sprime>1
           U_tplus1sprime=1; % U is constrained between 0 and 1
        end
        if StateMatrix(ActIn,2) == 1 
            Q_sa(gammaIdx,1,1,ActIn)= U_tplus1sprime-cost(ActIn)+gamma*Value(ActIn,2); % m is the maximum reward value
        else
            Q_sa(gammaIdx,1,1,ActIn)= -cost(ActIn)+gamma*Value(ActIn,2);
        end
    end
    Value(1,1)=max(Q_sa(gammaIdx,1,1,:));
    
    
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