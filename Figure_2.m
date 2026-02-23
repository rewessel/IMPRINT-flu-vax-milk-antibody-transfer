%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter to samples from  2 weeks postpartum & mom was vaccinated during
% pregnancy compared to no vaccine at all
samples = data(data.Time_point==1 & (data.Vax_timing_bin==1 |...
    data.Vax_timing_bin==3 | data.Vax_timing_bin==4 |...
    data.Vax_timing_bin==5),:); 

myColors = [245, 245, 245;166, 97, 26]/255;

X = log10(table2array(samples(:,6:284)));
Y = table2array(samples(:,5));
Y1 = zeros(height(Y),2);
Y1(Y==3 | Y==4 | Y==5,1) = 1; % first column is 1 if vaccinated at any point in pregnancy
Y1(Y==1,2) = 1; % second column is 1 if not vaccinated during pregnancy

Fig2 = PLSDA_main(X,Y1,2,myVarNames,{'yes',1,0.8},'orthogonal','',{'kfold',10},1000,{'vax','no vax'},myColors)

[~,Fig2.CV_pval]=ttest2(Fig2.CV_accuracy,Fig2.CV_accuracy_perm);