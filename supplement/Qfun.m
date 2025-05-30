% for HyperbolicDiscounting simulations
function Q_sa = Qfun(T,k,r_0,p,q,m,n,gammaVec,Deltas,StateVector,ActVector,StateMatrix)
%Deltas = 0.1;
%StateVector =  (0:Deltas:1)';
%StateMatrix =  nan(length(StateVector),T+1);
%StateMatrix(1,1)=0;
%StateMatrix(:,2:end)=repmat(StateVector,1,T);
%ActVector = StateVector/k; 
cost = p*ActVector.^q; 

Q_sa = nan(length(gammaVec),T,length(StateVector),length(ActVector));

for gammaIdx = 1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
%% First calculate Q(s,a) at T
% because the reward is delayed and only given at T+1, better to separate t<T and T.  
% If the reward is instantaneous, then we move this calculation to 'for loop'.
Value       =nan(length(StateVector),T+1);
Value(:,T+1)=0;

for stIn = 1:size(StateMatrix,1)
    NextStateVector = StateMatrix(stIn,T)+round(k*ActVector/Deltas)*Deltas;
    AvailActIndex = find(NextStateVector<=1);
    % for fun ActVector(1)=0;
    U_Tplus1 = m*StateMatrix(stIn,T)^n; % reward
    if U_Tplus1>1
        U_Tplus1=1; % U is constrained between 0 and 1
    end
    Q_sa(gammaIdx,T,stIn,1)= r_0+U_Tplus1+gamma*Value(stIn,T+1);
    % for work ActVector(>1)
    for ActIn = 2: length(AvailActIndex)
        U_Tplus1 = m*NextStateVector(ActIn)^n; % reward
        if U_Tplus1>1
            U_Tplus1=1; % U is constrained between 0 and 1
        end
        Q_sa(gammaIdx,T,stIn,ActIn)= -cost(ActIn)+U_Tplus1+gamma*Value(abs(NextStateVector(ActIn)-StateVector)<0.00001,T+1);
    end 
    Value(stIn,T)=nanmax(Q_sa(gammaIdx,T,stIn,:));
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
        % for fun ActVector(1)=0;
        Q_sa(gammaIdx,t,stIn,1)= r_0+gamma*Value(stIn,t+1);
        % for work
        for ActIn = 2: length(AvailActIndex)
            Q_sa(gammaIdx,t,stIn,ActIn)= -cost(ActIn)+gamma*Value(abs(NextStateVector(ActIn)-StateVector)<0.00001,t+1);
        end
        Value(stIn,t)=max(Q_sa(gammaIdx,t,stIn,:));

    end
end

end

end