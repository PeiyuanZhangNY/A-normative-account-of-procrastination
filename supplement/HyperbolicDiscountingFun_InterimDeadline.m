function [OptActSequence,OptStateSequence,net_earning] = HyperbolicDiscountingFun_InterimDeadline(T,k,r_0,p,q,m,n,k_DR,N_milestone,penalty_unit)

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

reward_interval_list = cumsum(ones(1,N_milestone)*1/N_milestone); %[1/3,2/3,1]
interim_deadline_list = cumsum(ones(1,N_milestone)*round(T/N_milestone)); % [7,14,21]

%% calculate Q_sa, here, very similar with function OptActStateSeq_interimdeadline with a few modifications marked

cost = p*ActVector.^q; 
gammaVec = TransGammaVector;
Q_sa = nan(length(gammaVec),T,length(StateVector),length(ActVector));

for gammaIdx = 1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    
    %% First calculate Q(s,a) at T 
    Value       =nan(length(StateVector),T+1);
    Value(:,T+1)=0;
    
    for stIn = 1:size(StateMatrix,1)
    NextStateVector = StateMatrix(stIn,T)+round(k*ActVector/Deltas)*Deltas;
    AvailActIndex = find(NextStateVector<=1);
    %AvailNextState = NextStateVector(AvailActIndex);
    % for fun ActVector(1)=0; final state is StateMatrix(stIn,T) 
    U_Tplus1 = m*StateMatrix(stIn,T)^n;
    if U_Tplus1>1
        U_Tplus1=1; % U is constrained between 0 and 1
    end    
    % if the final state is <1/3, -0.03*3 (penalty for delayed), if the
    % final state is 1/3=< <2/3, -0.03*2, if the final state is 2/3=< <1,
    % -0.03, if the final state is 1. no penalty. 
    listminusstate = reward_interval_list - StateMatrix(stIn,T); 
    closestindexOnLeftOfState = sum(listminusstate<=0); %
    % if closestindexOnLeftOfState=0, final state <1/3, if 1, final state
    % is between 1/3 and 2/3, if 2, final state between 2/3 and 1, if 3,
    % final state=1.
    penalty_loss_fun = (N_milestone-closestindexOnLeftOfState)*penalty_unit;
    
    Q_sa(gammaIdx,T,stIn,1)= r_0+U_Tplus1-penalty_loss_fun+gamma*Value(stIn,T+1);
    % if work;
    % StateMatrix(stIn,T+1) is the result of reward function
    for ActIn = 2: length(AvailActIndex)
        U_Tplus1 = m*NextStateVector(ActIn)^n;
        if U_Tplus1>1
            U_Tplus1=1; % U is constrained between 0 and 1
        end
        listminusstate = reward_interval_list - NextStateVector(ActIn); 
        closestindexOnLeftOfState = sum(listminusstate<=0);
        penalty_loss_work = (N_milestone-closestindexOnLeftOfState)*penalty_unit;
        
        Q_sa(gammaIdx,T,stIn,ActIn)= -cost(ActIn)+U_Tplus1-penalty_loss_work+gamma*Value(stIn+ActIn-1,T+1);
        %NextStateVector(QIn) is the result of reward function
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
        
        %%%%%%%%%% for fun 
        % when t<1/3*T (7)
        if t<interim_deadline_list(1)
           penalty_loss_fun = 0;
           penalty_loss_work = 0;
        end
        % when t between 1/3*T (7) and 2/3*T (14)
        if t>=interim_deadline_list(1) &&  t<interim_deadline_list(2)
            % for fun
            if StateMatrix(stIn,t)<1/N_milestone 
                penalty_loss_fun = penalty_unit;
            else
                penalty_loss_fun=0;
            end
%             % for work
%             if NextStateVector(QIn)<1/N_milestone
%                 penalty_loss_work = penalty_unit;
%             else
%                 penalty_loss_work = 0;
%             end           
        end 
        % when t between 2/3*T (14) and T (21)
        if t>=interim_deadline_list(2) &&  t<interim_deadline_list(3)
            % for fun
            if StateMatrix(stIn,t)<1/N_milestone 
                penalty_loss_fun = 2*penalty_unit;
            elseif StateMatrix(stIn,t)<2/N_milestone
                penalty_loss_fun = penalty_unit;
            else
                penalty_loss_fun=0;
            end
%             % for work
%             if NextStateVector(QIn)<1/N_milestone 
%                 penalty_loss_work = 2*penalty_unit;
%             elseif NextStateVector(QIn)<2/N_milestone
%                 penalty_loss_work = penalty_unit;
%             else
%                 penalty_loss_work=0;
%             end           
        end
        Q_sa(gammaIdx,t,stIn,1)= r_0-penalty_loss_fun+gamma*Value(stIn,t+1);
        
        %%%%%%%%%%% for work
        for ActIn = 2: length(AvailActIndex)
        % when t<1/3*T (7)
        if t<interim_deadline_list(1)
           penalty_loss_fun = 0;
           penalty_loss_work = 0;
        end
     
        % when t between 1/3*T (7) and 2/3*T (14)
        if t>=interim_deadline_list(1) &&  t<interim_deadline_list(2)
        % for work
            if NextStateVector(ActIn)<1/N_milestone
                penalty_loss_work = penalty_unit;
            else
                penalty_loss_work = 0;
            end  
        end
        
        % when t between 2/3*T (14) and T (21)
        if t>=interim_deadline_list(2) &&  t<interim_deadline_list(3)
        % for work
            if NextStateVector(ActIn)<1/N_milestone 
                penalty_loss_work = 2*penalty_unit;
            elseif NextStateVector(ActIn)<2/N_milestone
                penalty_loss_work = penalty_unit;
            else
                penalty_loss_work=0;
            end    
        end

        Q_sa(gammaIdx,t,stIn,ActIn)= -cost(ActIn)-penalty_loss_work+gamma*Value(stIn+ActIn-1,t+1);
        end
        
        % if there are multiple maximums, select the first maximum value
        % and index
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

%% net earning
net_earning = 0; % the utility gained in the end - cumulative penalty_loss along the way
% when t between 1/3*T (7) and 2/3*T (14), if <1/3, then lose penalty_unit
for t = interim_deadline_list(1):interim_deadline_list(2)-1
    if OptActSequence(t)<1/N_milestone
        net_earning = net_earning -penalty_unit;
    end
end

% when t between 2/3*T (14) and T (21), if <1/3, then lose 2*penalty_unit,
% if <2/3 but >=1/3, lose penalty_unit
cumsum_OptActSequence = cumsum(OptActSequence);
for t = interim_deadline_list(2):T-1
    if cumsum_OptActSequence(t)<1/N_milestone
        net_earning = net_earning -2*penalty_unit;
    elseif cumsum_OptActSequence(t)<2/N_milestone
        net_earning = net_earning -penalty_unit;
    end
end
% when t=T
listminusstate = reward_interval_list - cumsum_OptActSequence(T); 
closestindexOnLeftOfState = sum(listminusstate<=0); %
net_earning = net_earning-(N_milestone-closestindexOnLeftOfState)*penalty_unit;
net_earning = net_earning+m*cumsum_OptActSequence(T)^n;

end