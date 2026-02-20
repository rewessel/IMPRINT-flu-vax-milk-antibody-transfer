%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter to samples where mom got vaccinated during pregnancy
samples = data(data.Delta_weeks_bin==1 & (data.trimester_at_vax_First_Trimester==1 |...
    data.trimester_at_vax_Second_Trimester==1 | data.trimester_at_vax_Third_Trimester==1),:); 

X = log10(table2array(samples(:,4:282)));
Y = table2array(samples(:,[285,287,288]));

myColors = [223, 194, 125;245, 245, 245;128, 205, 193]/255;

Fig4 = PLSDA_main(X,Y,3,myVarNames,{'yes',0.3,0.8},'orthogonal','',{'kfold',10},...
    1000,{'1st Tri','2nd Tri','3rd Tri'},myColors)

Ygroup = [ones(1000,1);2*ones(1000,1);3*ones(1000,1);4*ones(1000,1)];
[~,~,stats] = kruskalwallis([Fig4.CV_accuracy(:,1); ...
    Fig4.CV_accuracy(:,2); ...
    Fig4.CV_accuracy(:,3); ...
    Fig4.CV_accuracy_perm'],Ygroup,'off');

c=multcompare(stats,'CriticalValueType','dunn-sidak','Display','off');

Fig4.CV_pval=c(:,6);