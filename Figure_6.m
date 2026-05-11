%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLSDA - compare breastmilk antibodies in babies who did or did not get an infection in first 6 months
%% FIGURE S7C-D

ids_to_exclude = [119:601,621,628];
ids_to_exclude = append('0',string(ids_to_exclude));

samples = combined_data(~contains(combined_data.Subject_ID,ids_to_exclude),:);
samples(30,7:end) = mean(samples(30:31,7:end));
samples = samples([1:30,32:end],contains(combined_data.Properties.VariableNames,'milk'));

myVarNames = samples.Properties.VariableNames;
myVarNames = strrep(myVarNames,'_','/');

X = table2array(samples);
Y = combined_data.flu_exposure_1_age_days(~contains(combined_data.Subject_ID,ids_to_exclude));
Y = Y([1:30,32:end]);

Y(Y<=180)=1; Y(Y>180)=0;
Y = logical(Y); Y = double([Y ~Y]);

c = sum(X~=0,1);
mdlVarNames = myVarNames(c>height(X)*0.1);
X = X(:,c>height(X)*0.1);
% mdlVarNames = myVarNames;

[varNames,ia] = run_elastic_net(zscore(log10(X+1)), Y, mdlVarNames, 'minMSE', 0.1, 100, 0.6, 5);
X = X(:,ia);
mdlVarNames = mdlVarNames(ia);

Fig5_milk = PLSDA_main(log10(X+1),Y,2,mdlVarNames,{'no',0.1,0.6},'orthogonal','',{'kfold',5},100,{'Flu+','Flu-'},myColors)

%% PLSDA - compare serum antibodies in babies who did or did not get an infection in first 6 months
%% FIGURE S7C-D

ids_to_exclude = [119:601,621,628];
ids_to_exclude = append('0',string(ids_to_exclude));

samples = combined_data(~contains(combined_data.Subject_ID,ids_to_exclude),:);
samples(30,7:end) = mean(samples(30:31,7:end));
samples = samples([1:30,32:end],contains(combined_data.Properties.VariableNames,'serum') | contains(combined_data.Properties.VariableNames,'milk'));

myVarNames = samples.Properties.VariableNames;
myVarNames = strrep(myVarNames,'_','/');

X = table2array(samples);
Y = combined_data.flu_exposure_1_age_days(~contains(combined_data.Subject_ID,ids_to_exclude));
Y = Y([1:30,32:end]);

Y(Y<=180)=1; Y(Y>180)=0;
Y = logical(Y); Y = double([Y ~Y]);

% filter low variance features to make it a little easier for elastic net
% to select serum features.
% myVarNames = myVarNames(var(log10(X))>0.15);
% X = X(:,var(log10(X))>0.15);

c = sum(X~=0,1);
mdlVarNames = myVarNames(c>height(X)*0.1);
X = X(:,c>height(X)*0.1);

[varNames,ia] = run_elastic_net(zscore(log10(X+1)), Y, mdlVarNames, 'minMSE', 0.1, 100, 0.6, 5);

X = X(:,ia);
mdlVarNames = mdlVarNames(ia);

Fig5_serum = PLSDA_main(log10(X+1),Y,2,mdlVarNames,{'no',0.1,0.6},'orthogonal','',{'kfold',5},100,{'Flu+','Flu-'},myColors)

%% manually combine milk and serum features into a single model

X = [Fig5_serum.XpreZ Fig5_milk.XpreZ];
Y = Y;
myVarNames = [Fig5_serum.varNames Fig5_milk.varNames];
Fig5_combined = PLSDA_main(X,Y,2,myVarNames,{'no',[],[]},'orthogonal','',{'kfold',5},1000,{'Flu+','Flu-'},myColors)

%% compare ROC curves for serum and breast milk models

figure
plot(Fig5_serum.ROC(:,1),Fig5_serum.ROC(:,2),'lineWidth',2); hold on
plot(Fig5_milk.ROC(:,1),Fig5_milk.ROC(:,2),'lineWidth',2)
plot(Fig5_combined.ROC(:,1),Fig5_combined.ROC(:,2),'lineWidth',2)
plot([0 1],[0 1],'k--')

xlabel('False Positive Rate'); ylabel('True Positive Rate')
legend(append('Serum (AUC=',string(round(Fig5_serum.AUC,2)),')'),...
    append('Milk (AUC=',string(round(Fig5_milk.AUC,2)),')'),...
    append('Combined (AUC=',string(round(Fig5_combined.AUC,2)),')'))

%% Create LASSO network for combined model
n=1;
fulldata = log10(combined_data(~contains(combined_data.Subject_ID,ids_to_exclude),7:537)+1);
fulldata(30,:) = mean(fulldata(30:31,:));
fulldata = fulldata([1:30,32:end],:);

fulldata_names = fulldata.Properties.VariableNames;
[rho,pval]=corr(Fig5_combined.XpreZ,table2array(fulldata),'type','Pearson');

Fig5_combined.varNames = strrep(Fig5_combined.varNames,'/','_')

for i = 1:width(Fig5_combined.XpreZ)
    for j = 1:width(fulldata)

    source_table(n) = Fig5_combined.varNames(i);
    target_table(n) = fulldata_names(j);
    rho_table(n) = rho(i,j);
    pval_table(n) = pval(i,j);
    n=n+1;
    end
end
pval_corrected = pval_table*length(pval_table); %apply bonferroni correction
LASSO_network = table(source_table',target_table',rho_table',pval_table',pval_corrected','VariableNames',{'LASSO_Feature','Correlate','rho','pval','pval_Bonferroni'});

    % remove self-correlations
LASSO_network = LASSO_network(~strcmp(LASSO_network.LASSO_Feature,LASSO_network.Correlate),:);

