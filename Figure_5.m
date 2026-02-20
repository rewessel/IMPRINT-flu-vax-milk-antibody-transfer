%% PLSDA modeling - comparing the effect of breast milk and serum antibodies at predicting future flu infections
cd 'Z:/users/wesselr/data/IMPRINT/'
data = readtable('IMPRINT_combined_serum_breastmilk_11102025.xlsx');

addpath('C:\Users\remzi\OneDrive\Documents\GitHub\PLSR-DA_MATLAB')
matvax = readtable('G_ALTER_Matvax_20250326_REW.xlsx');
metadata = readtable('G_ALTER_META_20250326_REW.csv');

data.Properties.VariableNames = strrep(data.Properties.VariableNames,'_','-');

myColors = [88 110 245;185 122 245;203 193 245;118 89 245; 86 158 245;219 86 245]/255;

% this is case control matching to validate PLSDA models in Fig 6
metadata_short = metadata(metadata.AGE_AT_BLOOD_DAYS_USE==0|metadata.AGE_AT_BLOOD_DAYS_USE==1,:);

% add a step to exclude children with no opportunity for flu exposure
ids_to_exclude = [119:601,621,628];
ids_to_exclude = append('0',string(ids_to_exclude));
metadata_short = metadata_short(~contains(metadata_short.Subject_ID,ids_to_exclude),:);

% extract flu positive and flu negative IDs
flu_pos = metadata_short(metadata_short.flu_exposure_1_age_days<=180 & ~strcmp(metadata_short.Aliquot_Type,'serum_m'),:);
flu_neg = metadata_short(metadata_short.flu_exposure_1_age_days>180 & ~strcmp(metadata_short.Aliquot_Type,'serum_m'),:);

% subset data based on exposure opportunity
data = data(~contains(data.SubjectID,ids_to_exclude),:);

%% PLSDA - compare breastmilk antibodies in babies who did or did not get an infection in first 6 months
%% FIGURE S7C-D

samples = [data.SubjectID data(:,contains(data.Properties.VariableNames,'milk'))]; 

samples.Properties.VariableNames{1} = 'Subject_ID';
myVarNames = samples.Properties.VariableNames(2:end);

[~,ia,ib]=intersect(samples.Subject_ID,metadata_short.Subject_ID);

joindata = join(samples(ia,:),metadata_short(ib,[1,17]),'Keys','Subject_ID');

X = table2array(joindata(:,2:280));
Y = table2array(joindata(:,281));

Y(Y<=180)=1; Y(Y>180)=0;
Y = logical(Y); Y = double([Y ~Y]);

[varNames,ia] = run_elastic_net(zscore(log10(X)), Y, myVarNames, 'minMSE', 0.1, 500, 0.6, 5);

X = X(:,ia);
myVarNames = myVarNames(ia);

Fig5_milk = PLSDA_main(log10(X),Y,2,myVarNames,{'no',nan,nan},'orthogonal','',{'kfold',10},1000,{'Flu+','Flu-'},myColors)

%% PLSDA - compare serum antibodies in babies who did or did not get an infection in first 6 months
%% FIGURE S7C-D

samples = [data.SubjectID data(:,contains(data.Properties.VariableNames,'serum'))]; % filter to samples from a certain time point
samples.Properties.VariableNames{1} = 'Subject_ID';
myVarNames = samples.Properties.VariableNames(2:end);

[~,ia,ib]=intersect(samples.Subject_ID,metadata_short.Subject_ID);

joindata = join(samples(ia,:),metadata_short(ib,[1,17]),'Keys','Subject_ID');

X = table2array(joindata(:,2:253));
Y = table2array(joindata(:,254));

Y(Y<=180)=1; Y(Y>180)=0;
Y = logical(Y); Y = double([Y ~Y]);

% filter low variance features out to help with elastic net
myVarNames = myVarNames(var(log10(X))>0.15);
X = X(:,var(log10(X))>0.15);

[varNames,ia] = run_elastic_net(zscore(log10(X)), Y, myVarNames, 'minMSE', 0.05, 500, 0.6, 5);

X = X(:,ia);
myVarNames = myVarNames(ia);
Fig5_serum = PLSDA_main(log10(X),Y,2,myVarNames,{'no',nan,nan},'orthogonal','',{'kfold',10},1000,{'Flu+','Flu-'},myColors)

%% manually combine milk and serum features into a single model

X = [Fig5_serum.XpreZ Fig5_milk.XpreZ];
Y = Y;
myVarNames = [Fig5_serum.varNames Fig5_milk.varNames];
Fig5_combined = PLSDA_main(X,Y,2,myVarNames,{'no',[],[]},'orthogonal','',{'kfold',10},1000,{'Flu+','Flu-'},myColors)

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
fulldata = log10(data(:,2:532));
fulldata_names = fulldata.Properties.VariableNames;
[rho,pval]=corr(Fig5_combined.XpreZ,table2array(fulldata),'type','Pearson');
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

