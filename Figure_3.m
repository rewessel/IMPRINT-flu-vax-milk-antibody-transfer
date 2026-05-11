%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter to samples from  2 weeks postpartum & mom was vaccinated during
% pregnancy compared to no vaccine at all
samples = data(data.Time_point==1 & (data.Vax_timing_bin==1 |...
    data.Vax_timing_bin==3 | data.Vax_timing_bin==4 |...
    data.Vax_timing_bin==5),:); 
% samples = data(data.Time_point==1 & (data.Vax_timing_bin==1 |...
%     data.Vax_timing_bin==3),:); 
% samples = data(data.Time_point==1 & (data.Vax_timing_bin==1 |...
%     data.Vax_timing_bin==4 | data.Vax_timing_bin==5),:); 


samples(8,6:end) = mean(samples(8:9,6:end)); 
samples = samples([1:8,10:end],:); % remove duplicated sample

myColors = [245, 245, 245;166, 97, 26]/255;

X = log10(table2array(samples(:,6:284))+1);
Y = table2array(samples(:,5));
Y1 = zeros(height(Y),2);
Y1(Y==3 | Y==4 | Y==5,1) = 1; % first column is 1 if vaccinated at any point in pregnancy
Y1(Y==1,2) = 1; % second column is 1 if not vaccinated during pregnancy

% remove features that are more than 90% zeros

c = sum(X~=0,1);
mdlVarNames = myVarNames(c>height(X)*0.1);
X = X(:,c>height(X)*0.1);
% mdlVarNames = myVarNames;

% X = X([1:43,45:49],:);
% Y1 = Y1([1:43,45:49],:);

% [varNames,ia] = run_elastic_net(zscore(X), Y1, mdlVarNames, 'minMSE',0.1, 100, 0.8, 5);
Fig2 = PLSDA_main(X,Y1,2,mdlVarNames,{'yes',0.1,0.8},'orthogonal','',{'kfold',5},100,{'vax','no vax'},myColors)
% Fig2 = PLSDA_main(X,Y1,2,mdlVarNames,{'yes',1,0.8},'orthogonal','',{'kfold',10},1000,{'vax','no vax'},myColors)

[~,Fig2.CV_pval]=ttest2(Fig2.CV_accuracy,Fig2.CV_accuracy_perm);