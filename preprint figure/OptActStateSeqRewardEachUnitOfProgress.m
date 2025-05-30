function [ OptActSequence,OptStateSequence,Value ] = OptActStateSeqRewardEachUnitOfProgress( Deltas,T,k,r_0,p,q,gamma,m,n )
% same with the function with the original name OptActStateSeqNoDelay.m
% I just change the name to be aligned with the paper

% this is the version without delayed feedback U_t = alpha*s_t^beta  
% deltas = k*a
%Deltas      =  0.1;
StateVector =  (0:Deltas:1)';
%T           =  10;
StateMatrix =  nan(length(StateVector),T+1);
StateMatrix(1,1)=0;
StateMatrix(:,2:end)=repmat(StateVector,1,T);

% action deltas=k*a
%k = 0.5; 
ActVector = StateVector/k; % because state unit is fixed here, so k<=1,
% if k>1, by having an action ActVector, the state will end up somewhere
% within the unit. 

% reward values
%r_0 = 0.1;
% R = s, could be a separate function, but for simplicity here, R=s. 
%p=2; q=6;
cost = p*ActVector.^q; 
% cost =  p*ActVector; 
% cost(q:end)=2;
% myopia parameter
%gamma = 0.9;


Value       =nan(length(StateVector),T+1);
Value(:,T+1)=0;
OptActIn    =nan(length(StateVector),T); % choose the action with index in ActVector.

for t = T:-1:1
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
        U_tplus1s = m*StateMatrix(stIn,t+1)^n;
        if U_tplus1s>1
            U_tplus1s=1; % U is constrained between 0 and 1
        end
        U_ts = m*StateMatrix(stIn,t)^n;
        if U_ts>1
            U_ts=1; % U is constrained between 0 and 1
        end
        % mathematically U_tplus1s always = U_ts
        Q(1)= r_0+U_tplus1s-U_ts+gamma*Value(stIn,t+1); 
        % for work
        
        for QIn = 2: size(Q,1)
            U_tplus1sprime = m*StateMatrix(stIn+QIn-1,t+1)^n;
            if U_tplus1sprime>1
                U_tplus1sprime=1; % U is constrained between 0 and 1
            end
            
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

