%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter to samples where mom got vaccinated during pregnancy
samples = data(data.Time_point==1 & (data.Vax_timing_bin==3 |...
    data.Vax_timing_bin==4 | data.Vax_timing_bin==5),:); 

X = log10(table2array(samples(:,6:284)));
Y = categorical(samples.Vax_timing_bin);
Y = onehotencode(Y,2);
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