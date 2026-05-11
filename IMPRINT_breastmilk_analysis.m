%% IMPRINT breast milk manuscript analysis & multivariate modeling
cd 'C:\Users\remzi\OneDrive\Documents\GitHub\Flu-vaccination-breastmilk-transfer'
addpath('C:\Users\remzi\OneDrive\Documents\GitHub\PLSR-DA_MATLAB')

data = readtable('Milk_serology_v3.xlsx');
% data = readtable('Milk_serology.xlsx');

data = data(:,~contains(data.Properties.VariableNames,'Ebola') &...
    ~contains(data.Properties.VariableNames,'RSV') & ...
    ~contains(data.Properties.VariableNames,'EBV'));
matvax = readtable('Maternal_flu_vaccine_records.xlsx');

% data subset to combined serum and milk features, with only features that
% had both cord serum and milk @ 2 weeks postpartum
% this dataframe has also excluded kids who did not have an opportunity to
% be exposed to flu (i.e., their first 6 months of life did not coincide
% with that year's flu season).
combined_data = readtable('Fig6_combined_serum_milk_v2.xlsx');

myVarNames = data.Properties.VariableNames(6:284);
myVarNames = regexprep(myVarNames,'_','-');

myColors = [88 110 245;185 122 245;203 193 245;118 89 245; 86 158 245;219 86 245]/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Figure_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 - CROSS-VALIDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = sum(X~=0,1);
mdlVarNames = myVarNames(c>height(X)*0.1);
X = X(:,c>height(X)*0.1);

% [Fig2_CV.cv_true, Fig2_CV.cv_perm, Fig2_CV.cv_rand, ~] = ...
%     crossValidate(X,Y1,myVarNames,10,2,5,1,0.8,'not-multilevel',10)
[Fig2_CV.cv_true, Fig2_CV.cv_perm, Fig2_CV.cv_rand, ~,...
    Fig2_CV.Ytrue_predicted, Fig2_CV.Yperm_predicted, Fig2_CV.Yrand_predicted] = ...
    crossValidate(X,Y1,mdlVarNames,height(X),2,50,0.1,0.8,'not-multilevel',100)
% 04/30 this generated one good run - worth doing a full set of runs with
% more lasso tests?
figure;
swarmchart(ones(1,length(Fig2_CV.cv_true)),Fig2_CV.cv_true); hold on
swarmchart(2*ones(1,length(Fig2_CV.cv_rand)),Fig2_CV.cv_rand)
swarmchart(3*ones(1,length(Fig2_CV.cv_perm)),Fig2_CV.cv_perm)
ylabel('Accuracy');xticks([1:3]);xticklabels({'true','rand','perm'})
% If this fails, try the same protocol with the EBOLA normalized data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Figure_3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Figure_4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4 - CROSS-VALIDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Fig4_CV.cv_true, Fig4_CV.cv_perm, Fig4_CV.cv_rand, ~] = ...
    crossValidate(X,Y,myVarNames,5,3,10,1,0.8,'not-multilevel',10)

figure;
swarmchart(ones(1,length(Fig4_CV.cv_true)),Fig4_CV.cv_true); hold on
swarmchart(2*ones(1,length(Fig4_CV.cv_rand)),Fig4_CV.cv_rand)
swarmchart(3*ones(1,length(Fig4_CV.cv_perm)),Fig4_CV.cv_perm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Figure_5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5 - CROSS-VALIDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Fig5_CV.cv_true, Fig5_CV.cv_perm, Fig5_CV.cv_rand, ~] = ...
    crossValidate(X,Y,mdlVarNames,5,2,10,0.1,0.6,'not-multilevel',100)