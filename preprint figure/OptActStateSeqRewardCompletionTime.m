function [ OptActSequence,OptStateSequence,Value ] = OptActStateSeqRewardCompletionTime( Deltas,T,k,r_0,p,q,gamma,m,n )
% reward upon task completion

% deltas = k*a
% Deltas      =  0.01;
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
% m=1;n=1;
% gamma = 1;

%% First calculate Q(s,a) at T 
Value       =nan(length(StateVector),T+1);
Value(:,T+1)=0;
OptActIn    =nan(length(StateVector),T); % choose the action with index in ActVector.

% because the reward given at T+1 is different from that given before T+1.

% calculate the current state<1 (not the final state)
for stIn = 1:size(StateMatrix,1)-1 % The last state is 1, then it should have been rewarded at T-1, so should be separately analyzed
    NextStateVector = StateMatrix(stIn,T)+k*ActVector;
    AvailActIndex = find(NextStateVector<=1);
    Q = nan(length(AvailActIndex),1);
    
    % for fun ActVector(1)=0;
    U_Tplus1 = m*StateMatrix(stIn,T+1)^n;
    if U_Tplus1>1
        U_Tplus1=1; % U is constrained between 0 and 1
    end
    Q(1)= r_0+U_Tplus1+gamma*Value(stIn,T+1);
    
    % for work
    for QIn = 2: size(Q,1)
        U_tplus1sprime = m*StateMatrix(stIn+QIn-1,T+1)^n;
        if U_tplus1sprime>1
            U_tplus1sprime=1; % U is constrained between 0 and 1
        end
        Q(QIn)= -cost(QIn)+U_tplus1sprime+gamma*Value(stIn+QIn-1,T+1);

    end 
    [Value(stIn,T),OptActIn(stIn,T)]=max(Q);
end
% if the current state is 1 
OptActIn(size(StateMatrix,1),T) = 1; % choose for fun (only choice when you already completed the task at T)
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
        Q = nan(length(AvailActIndex),1);
        % for fun ActVector(1)=0;
        Q(1)= r_0+gamma*Value(stIn,t+1); 
        
        % for work        
        for QIn = 2: size(Q,1)
            U_tplus1sprime = m*StateMatrix(stIn+QIn-1,t+1)^n;
            if U_tplus1sprime>1
               U_tplus1sprime=1; % U is constrained between 0 and 1
            end
            if StateMatrix(stIn+QIn-1,t+1) == 1 
                Q(QIn)= U_tplus1sprime-cost(QIn)+gamma*Value(stIn+QIn-1,t+1); % m is the maximum reward value
            else
                Q(QIn)= -cost(QIn)+gamma*Value(stIn+QIn-1,t+1);
            end
        end
        [Value(stIn,t),OptActIn(stIn,t)]=max(Q);
        %Value(stIn,t)=max(Q);
        %tempOptActIn = find(abs(Q - Value(stIn,t))<eps);
        %tempSeq = randperm(length(tempOptActIn));
        %OptActIn(stIn,t) = tempSeq(1);
    end
    
    % if the current state is 1 

    OptActIn(size(StateMatrix,1),t) = 1; % choose for fun (only choice when you already completed the task at T)
    Value(size(StateMatrix,1),t) = r_0+gamma*Value(size(StateMatrix,1),t+1);
    
end

%% last calculate Q(s,a) at t=1
NextStateVector = StateMatrix(1,1)+k*ActVector;
AvailActIndex = find(NextStateVector<=1);
Q = nan(length(AvailActIndex),1);
% for fun ActVector(1)=0;
Q(1)= r_0+gamma*Value(1,2); 

% for work        
for QIn = 2: size(Q,1)
    U_tplus1sprime = m*StateMatrix(QIn,2)^n;
    if U_tplus1sprime>1
       U_tplus1sprime=1; % U is constrained between 0 and 1
    end
    if StateMatrix(QIn,2) == 1 
        Q(QIn)= U_tplus1sprime-cost(QIn)+gamma*Value(QIn,2); % m is the maximum reward value
    else
        Q(QIn)= -cost(QIn)+gamma*Value(QIn,2);
    end
end
[Value(1,1),OptActIn(1,1)]=max(Q);
  


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

