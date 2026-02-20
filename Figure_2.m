%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter to samples from  2 weeks postpartum & mom was vaccinated during
% pregnancy compared to no vaccine at all
samples = data(data.Delta_weeks_bin==1 & (data.trimester_at_vax_First_Trimester==1 |...
    data.trimester_at_vax_Second_Trimester==1 | data.trimester_at_vax_Third_Trimester==1 |...
    data.trimester_at_vax_No_Vaccine==1),:); 

myColors = [245, 245, 245;166, 97, 26]/255;

X = log10(table2array(samples(:,4:282)));
Y = table2array(samples(:,283:288));
Y1 = zeros(height(Y),2);
Y1(Y(:,3)==1 | Y(:,5)==1 | Y(:,6)==1,1) = 1; % first column is 1 if vaccinated at any point in pregnancy
Y1(Y(:,4)==1,2) = 1; % second column is 1 if not vaccinated during pregnancy

Fig2 = PLSDA_main(X,Y1,2,myVarNames,{'yes',1,0.8},'orthogonal','',{'kfold',10},1000,{'vax','no vax'},myColors)

[~,Fig2.CV_pval]=ttest2(Fig2.CV_accuracy,Fig2.CV_accuracy_perm);