function OptActSequence = HyperbolicDiscountingFun_RewardEachUnitOfProgress(T,k,r_0,p,q,m,n,k_DR)

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

%% calculate Q_sa, here, very similar with function OptActStateSeqRewardEachUnitOfProgress with a few modifications marked

cost = p*ActVector.^q; 
gammaVec = TransGammaVector;
Q_sa = nan(length(gammaVec),T,length(StateVector),length(ActVector));

for gammaIdx = 1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    
    Value       =nan(length(StateVector),T+1);
    Value(:,T+1)=0;
    
    for t = T:-1:1
    if t==1
        stInend = 1;
    else 
        stInend = size(StateMatrix,1);
    end
    for stIn = 1:stInend
        NextStateVector = StateMatrix(stIn,t)+k*ActVector;
        AvailActIndex = find(NextStateVector<=1);
        % for fun ActVector(1)=0;
        U_tplus1s = m*StateMatrix(stIn,t+1)^n;
        if U_tplus1s>1
            U_tplus1s=1; % U is constrained between 0 and 1
        end
        U_ts = m*StateMatrix(stIn,t)^n;
        if U_ts>1
            U_ts=1; % U is constrained between 0 and 1
        end
        % mathematically U_tplus1s always = U_ts
        Q_sa(gammaIdx,T,stIn,1)= r_0+U_tplus1s-U_ts+gamma*Value(stIn,t+1); 
        % for work
        
        for ActIn = 2: length(AvailActIndex)
            U_tplus1sprime = m*StateMatrix(stIn+ActIn-1,t+1)^n;
            if U_tplus1sprime>1
                U_tplus1sprime=1; % U is constrained between 0 and 1
            end
            
            Q_sa(gammaIdx,t,stIn,ActIn)= U_tplus1sprime-U_ts-cost(ActIn)+gamma*Value(stIn+ActIn-1,t+1);
        end
        Value(stIn,t)=nanmax(Q_sa(gammaIdx,t,stIn,:));
        %Value(stIn,t)=max(Q);
        %tempOptActIn = find(abs(Q - Value(stIn,t))<eps);
        %tempSeq = randperm(length(tempOptActIn));
        %OptActIn(stIn,t) = tempSeq(1);
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