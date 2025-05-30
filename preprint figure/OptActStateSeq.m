function [ OptActSequence,OptStateSequence,Value ] = OptActStateSeq( Deltas,T,k,r_0,p,q,gamma,m,n )
% this is the function that works for k<=1, this function resolves
% the issue that deltas and deltasa are different when k<1. 

StateVector =  (0:Deltas:1)';
StateMatrix =  nan(length(StateVector),T+1);
StateMatrix(1,1)=0;
StateMatrix(:,2:end)=repmat(StateVector,1,T);
ActVector = StateVector/k; 
cost = p*ActVector.^q; 

%% First calculate Q(s,a) at T 
% because the reward is only given at T+1, better to separate t<T and T.  
% If the reward is delievered along the task duration,
% then we put this calculation input 'for loop'.
Value       =nan(length(StateVector),T+1);
Value(:,T+1)=0;
OptActIn    =nan(length(StateVector),T); % choose the action with index in ActVector.

for stIn = 1:size(StateMatrix,1)
    NextStateVector = StateMatrix(stIn,T)+round(k*ActVector/Deltas)*Deltas;
    AvailActIndex = find(NextStateVector<=1);
    %AvailNextState = NextStateVector(AvailActIndex);
    Q = nan(length(AvailActIndex),1);
    % for fun ActVector(1)=0;
    U_Tplus1 = m*StateMatrix(stIn,T)^n;
    if U_Tplus1>1
        U_Tplus1=1; % U is constrained between 0 and 1
    end
    Q(1)= r_0+U_Tplus1+gamma*Value(stIn,T+1);
    % StateMatrix(stIn,T+1) is the result of reward function
    for QIn = 2: size(Q,1)
        U_Tplus1 = m*NextStateVector(QIn)^n;
        if U_Tplus1>1
            U_Tplus1=1; % U is constrained between 0 and 1
        end
        Q(QIn)= -cost(QIn)+U_Tplus1+gamma*Value(abs(NextStateVector(QIn)-StateVector)<0.00001,T+1);
        %NextStateVector(QIn) is the result of reward function
    end 
    [Value(stIn,T),OptActIn(stIn,T)]=max(Q);
%     % if there are multiple maximums, select the last maximum value and
%         % index
%         Value(stIn,T)=max(Q);
%         multipleidx = find(Q==Value(stIn,T));
%         OptActIn(stIn,T)=multipleidx(end);
end

%% Second, calculate Q(s,a) at t<T.

for t = T-1:-1:1
    if t==1
        stInend = 1;
    else 
        stInend = size(StateMatrix,1);
    end
    for stIn = 1:stInend
        NextStateVector = StateMatrix(stIn,t)+round(k*ActVector/Deltas)*Deltas;
        AvailActIndex = find(NextStateVector<=1);
        Q = nan(length(AvailActIndex),1);
        % for fun ActVector(1)=0;
        Q(1)= r_0+gamma*Value(stIn,t+1);
        % for work
        for QIn = 2: size(Q,1)
            Q(QIn)= -cost(QIn)+gamma*Value(abs(NextStateVector(QIn)-StateVector)<0.00001,t+1);
        end
        % if there are multiple maximums, select the first maximum value
        % and index
        [Value(stIn,t),OptActIn(stIn,t)]=max(Q); 
        
%         % if there are multiple maximums, select the last maximum value and
%         % index
%         Value(stIn,t)=max(Q);
%         multipleidx = find(Q==Value(stIn,t));
%         OptActIn(stIn,t)=multipleidx(end);

    end
end

%% find the optimal action sequence and its corresponding state sequence
% with state action trainsition rule 
OptStateSequence = nan(1,T+1);
OptStateSequence(1)=0;
OptStateSequence(2) = round(ActVector(OptActIn(1,1))*k/Deltas)*Deltas; 

OptActSeqenceIndex=nan(1,T);
OptActSeqenceIndex(1)=OptActIn(1,1);
for t = 2: T
    indexTemp=find(abs(StateMatrix(:,t)-OptStateSequence(t))<0.000001);
    %disp(['optimal actions at',num2str(t),'and',num2str(indexTemp)])
    
    OptActSeqenceIndex(t)=OptActIn(indexTemp,t);
    delta_s= round(ActVector(OptActSeqenceIndex(t))*k/Deltas)*Deltas;
    OptStateSequence(t+1)=OptStateSequence(t)+delta_s;
end
OptActSequence = ActVector(OptActSeqenceIndex); 


end

