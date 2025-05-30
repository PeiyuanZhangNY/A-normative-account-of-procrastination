function OptActSequence = HyperbolicDiscountingFun(T,k,r_0,p,q,m,n,k_DR)
%DeltaGammaVector = 0.001;
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


Q_saGamma = Qfun(T,k,r_0,p,q,m,n,TransGammaVector,Deltas,StateVector,ActVector,StateMatrix);

%% calculate Q_saK
Q_saK = nan(T,length(StateVector),length(ActVector));
OptActIn = nan(length(StateVector),T);

for t = 1:T
    for stIn = 1:size(StateMatrix,1)
        NextStateVector = StateMatrix(stIn,T)+round(k*ActVector/Deltas)*Deltas;
        AvailActIndex = find(NextStateVector<=1);
        
        for ActIn = 1: length(AvailActIndex)
            Q_saK(t,stIn,ActIn) = nansum(Q_saGamma(:,t,stIn,ActIn))*DeltaGammaVector;
        end
        [~,OptActIn(stIn,t)] = nanmax(Q_saK(t,stIn,:));
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