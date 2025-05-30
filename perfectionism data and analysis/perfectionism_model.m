% the data is perfectionism_rawdata sheet: remove missing data, columns: 0%/nothing	20%	40%	60%	80%	100%/all

% import data and load data
load('satisfactionrawdata.mat')

% convert satisfactionrawdata to percent e.g., 80 to 0.8
% fit satisfaction level with function m*x^a, if the fitting can not be
% done because of poor relationship, then mark the exponent vector NaN.
% a=1 neutral, a>1 perfectionist, a<1 not perfectionist. So log(a)=0 natural, log(a)>0 perfectionist, log(a)<0 not perfectionist.
SatisfactionLevel = satisfactionrawdata/100;
ft = fittype('m*x^a');
effortvec = [0,0.2,0.4,0.6,0.8,1]';
amvec = nan(size(SatisfactionLevel,1),2);  %[a,m]
for ii = 1: size(SatisfactionLevel,1)
    satisfvec = SatisfactionLevel(ii,:)';
    try 
        f=fit(effortvec,satisfvec,ft,'StartPoint',[1,1]); 
        amvec(ii,:) = coeffvalues(f); % [a,m]
    catch
        fprintf('loop number %d failed\n',ii)
    end    
end

exponentvec = amvec(:,1);
mvec = amvec(:,2);
logexponent = log(exponentvec);
histogram(logexponent)
nanmean(logexponent) % 0.4859 so average exponent is 1.6256 (exp(0.4859)), which means the average people tend to pursue perfectionism


% show the cases where exponent fitting failed: 
failedExponetfit = SatisfactionLevel(isnan(logexponent),:);
length(failedExponetfit) % 11
% these are the cases where the satisfaction level doesn't follow the
% monotonic increase with the increased final proportion completed. 

% predicted satisfaction level 
pred_satis_effort1 = 0;
pred_satis_effort2 = mvec.*effortvec(2).^exponentvec*100;
pred_satis_effort3 = mvec.*effortvec(3).^exponentvec*100;
pred_satis_effort4 = mvec.*effortvec(4).^exponentvec*100;
pred_satis_effort5 = mvec.*effortvec(5).^exponentvec*100;
pred_satis_effort6 = mvec*100;

% save the result 
save('perfectionsim_model.mat')