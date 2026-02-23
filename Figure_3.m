%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3C, Figure S6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear fc
data = sortrows(data,{'SubjectID','Delta_weeks'});
studyID = unique(data.SubjectID);

fc.IgG1 = nan(length(studyID),7);
fc.IgG2 = nan(length(studyID),7);
fc.IgG3 = nan(length(studyID),7);
fc.IgG4 = nan(length(studyID),7);
fc.IgM = nan(length(studyID),7);
fc.IgA1 = nan(length(studyID),7);
fc.IgA2 = nan(length(studyID),7);
fc.FcgR2A = nan(length(studyID),7);
fc.FcgR2B = nan(length(studyID),7);
fc.FcgR3A = nan(length(studyID),7);
fc.FcgR3B = nan(length(studyID),7);
fc.FcaR = nan(length(studyID),7);
fc.ADCD = nan(length(studyID),7);
fc.ADNP = nan(length(studyID),7);
fc.ADCP = nan(length(studyID),7);

% uncomment whichever antigen you wish to plot.
myAg = 'Kansas';
% myAg = 'Victoria';
% myAg = 'Phuket';
% myAg = 'Austria';


for i = 1:length(studyID)
    clear fc_temp fc_array
    current_patient = studyID(i);
    current_data = data(strcmp(data.SubjectID,current_patient),4:282);
    time_data = data.Delta_weeks(strcmp(data.SubjectID,current_patient));

    for j = 1:height(current_data)
        fc_temp(j,:) = current_data(j,:);
        time_points(j,i) = time_data(j);
    end

    varNames = fc_temp.Properties.VariableNames;
    fc_array = table2array(fc_temp);
    fc.IgG1(i,1:j) = mean(fc_array(:,contains(varNames,'IgG1') & contains(varNames,myAg)),2)';
    fc.IgG2(i,1:j) = mean(fc_array(:,contains(varNames,'IgG2') & contains(varNames,myAg)),2)';
    fc.IgG3(i,1:j) = mean(fc_array(:,contains(varNames,'IgG3') & contains(varNames,myAg)),2)';
    fc.IgG4(i,1:j) = mean(fc_array(:,contains(varNames,'IgG4') & contains(varNames,myAg)),2)';
    fc.IgM(i,1:j) = mean(fc_array(:,contains(varNames,'IgM') & contains(varNames,myAg)),2)';
    fc.IgA1(i,1:j) = mean(fc_array(:,contains(varNames,'IgA1') & contains(varNames,myAg)),2)';
    fc.IgA2(i,1:j) = mean(fc_array(:,contains(varNames,'IgA2') & contains(varNames,myAg)),2)';
    fc.FcgR2A(i,1:j) = mean(fc_array(:,contains(varNames,'FcgR2A') & contains(varNames,myAg)),2)';
    fc.FcgR2B(i,1:j) = mean(fc_array(:,contains(varNames,'FcgR2B') & contains(varNames,myAg)),2)';
    fc.FcgR3A(i,1:j) = mean(fc_array(:,contains(varNames,'FcgR3A') & contains(varNames,myAg)),2)';
    fc.FcgR3B(i,1:j) = mean(fc_array(:,contains(varNames,'FcgR3B') & contains(varNames,myAg)),2)';
    fc.FcaR(i,1:j) = mean(fc_array(:,contains(varNames,'FcaR') & contains(varNames,myAg)),2)';
    fc.ADCD(i,1:j) = mean(fc_array(:,contains(varNames,'ADCD') & contains(varNames,myAg)),2)';
    fc.ADNP(i,1:j) = mean(fc_array(:,contains(varNames,'ADNP') & contains(varNames,myAg)),2)';
    fc.ADCP(i,1:j) = mean(fc_array(:,contains(varNames,'ADCP') & contains(varNames,myAg)),2)';

end

time_points(time_points==0) = nan;

%% Plot the results
figure
fns = fieldnames(fc);
for i = 1:15
        nexttile; plotdata = log10(fc.(fns{i})');

    for j = 1:length(studyID)
            plot(time_points(1:3,j),plotdata(1:3,j),'.-','color',[0.1 0.1 0.1 0.2]); hold on
    end

    title(fns(i)); xlabel('Weeks postpartum'); ylabel('Log10(MFI)');xlim([0 52]) 

end

