function [ OptActSequence,OptStateSequence,net_earning ] = OptActStateSeq_interimdeadline( Deltas,T,k,r_0,p,q,gamma,m,n,N_milestone,penalty_unit )
% delayed reward with intermediate deadline: if not met, loose point immediately, but final reward is delayed.
% Total time=21 days. check every 1/3*21=7 days. linear reward
% when t=?*T, cumulative progress (proportion completed so far)=?, if met, nothing happened (just -cost), if not met, -cost-1% (Here I assume loose 0.01 grade immediately. People know if they are behind the deadline, they will loose grade 0.01 immediately); then when t=?*T+1, if met, nothing happened (-cost), if not met, -cost-1%; etc.  
%when t=?*T, continue to test if cumulative progress =?, if not pass, apply both penalty for the first deadline and second deadline, -2%, if ? is met, but ? is not met,  -1%, if ? is met, nothing happened. 
%When t=T, if cumulative progress<?, linear reward cumulative progress-2%, if ? is met, but ? is not met,  cumulative progress-1%, if 1 is met, reward=1. Here, it is reasonable to assume a linear reward function in this task, more progress is made, higher quality of the paper, higher grade (maximum grade is 1, 100%). 


StateVector =  (0:Deltas:1)';
StateMatrix =  nan(length(StateVector),T+1);
StateMatrix(1,1)=0;
StateMatrix(:,2:end)=repmat(StateVector,1,T);
ActVector = StateVector/k; 
cost = p*ActVector.^q; 

%N_milestone =3; % three texts for error correction 
reward_interval_list = cumsum(ones(1,N_milestone)*1/N_milestone); %[1/3,2/3,1]
interim_deadline_list = cumsum(ones(1,N_milestone)*round(T/N_milestone)); % [7,14,21]

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
    
    Q(1)= r_0+U_Tplus1-penalty_loss_fun+gamma*Value(stIn,T+1);
    % if work;
    % StateMatrix(stIn,T+1) is the result of reward function
    for QIn = 2: size(Q,1)
        U_Tplus1 = m*NextStateVector(QIn)^n;
        if U_Tplus1>1
            U_Tplus1=1; % U is constrained between 0 and 1
        end
        listminusstate = reward_interval_list - NextStateVector(QIn); 
        closestindexOnLeftOfState = sum(listminusstate<=0);
        penalty_loss_work = (N_milestone-closestindexOnLeftOfState)*penalty_unit;
        
        Q(QIn)= -cost(QIn)+U_Tplus1-penalty_loss_work+gamma*Value(abs(NextStateVector(QIn)-StateVector)<0.00001,T+1);
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
        Q(1)= r_0-penalty_loss_fun+gamma*Value(stIn,t+1);
        
        %%%%%%%%%%% for work
        for QIn = 2: size(Q,1)
        % when t<1/3*T (7)
        if t<interim_deadline_list(1)
           penalty_loss_fun = 0;
           penalty_loss_work = 0;
        end
     
        % when t between 1/3*T (7) and 2/3*T (14)
        if t>=interim_deadline_list(1) &&  t<interim_deadline_list(2)
        % for work
            if NextStateVector(QIn)<1/N_milestone
                penalty_loss_work = penalty_unit;
            else
                penalty_loss_work = 0;
            end  
        end
        
        % when t between 2/3*T (14) and T (21)
        if t>=interim_deadline_list(2) &&  t<interim_deadline_list(3)
        % for work
            if NextStateVector(QIn)<1/N_milestone 
                penalty_loss_work = 2*penalty_unit;
            elseif NextStateVector(QIn)<2/N_milestone
                penalty_loss_work = penalty_unit;
            else
                penalty_loss_work=0;
            end    
        end

        Q(QIn)= -cost(QIn)-penalty_loss_work+gamma*Value(abs(NextStateVector(QIn)-StateVector)<0.00001,t+1);
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

