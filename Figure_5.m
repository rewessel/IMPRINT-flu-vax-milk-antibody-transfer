%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter to samples where mom got vaccinated during pregnancy
% samples = data(data.Time_point==1 & (data.Vax_timing_bin==1 | data.Vax_timing_bin==3 |...
%     data.Vax_timing_bin==4 | data.Vax_timing_bin==5),:); 
samples = data(data.Time_point==1 & (data.Vax_timing_bin==3 |...
    data.Vax_timing_bin==4 | data.Vax_timing_bin==5),:); 

samples(8,6:end) = mean(samples(8:9,6:end)); 
samples = samples([1:8,10:end],:); % remove duplicated sample

X = log10(table2array(samples(:,6:284))+1);
Y = categorical(samples.Vax_timing_bin);
Y = onehotencode(Y,2);
myColors = [223, 194, 125;245, 245, 245;128, 205, 193;0 0 0]/255;

% remove features more than 25% zeroes
c = sum(X~=0,1);

mdlVarNames = myVarNames(c>height(X)*0.1);
X = X(:,c>height(X)*0.1);
% mdlVarNames = myVarNames;

% [varNames,ia] = run_elastic_net(zscore(X), Y, myVarNames, 'minMSE',0.6, 100, 0.8, 10);

Fig4 = PLSDA_main(X,Y,3,mdlVarNames,{'yes',0.6,0.8},'orthogonal','',{'kfold',5},...
    100,{'no vax','1st Tri','2nd Tri','3rd Tri'},myColors)


Ygroup = [ones(100,1);2*ones(100,1);3*ones(100,1);4*ones(100,1)];
[~,~,stats] = kruskalwallis([Fig4.CV_accuracy(:,1); ...
    Fig4.CV_accuracy(:,2); ...
    Fig4.CV_accuracy(:,3); ...
    Fig4.CV_accuracy_perm'],Ygroup,'off');

c=multcompare(stats,'CriticalValueType','dunn-sidak','Display','off');

Fig4.CV_pval=c(:,6);