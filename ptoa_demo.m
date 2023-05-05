%#######################################################################
%
%                    * PTOA DEMOgraphics Program *
%
%          M-File which reads the PTOA demographics data from the CSV 
%     file in the Results\ subdirectory:
%     PTOAStudyTheAdaptive-DataPull3May2023_DATA_2023-05-03_1359.csv.
%
%          The PTOA T1rho and T2* results are read and the demographics
%     are incorporated into the results.  The expanded results are
%     written to the MS-Excel spreadsheet:  mri_fitps_stat.xlsx in the
%     PTOA Results\mri_fitps\ subdirectory.
%
%     NOTES:  None.
%
%     05-May-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Read Demographics Data for All PTOA Subjects from CSV File
%
data = readtable(fullfile('Results', ...
     'PTOAStudyTheAdaptive-DataPull3May2023_DATA_2023-05-03_1359.csv'));
subjID = data.ptoa_subj;
subjID = strrep(subjID,'PTOA','');
subjID = strrep(subjID,'-1','');
%
% Get PTOA Subject Demographics
%
subjn = str2num(char(subjID));
sex = 2-(data(:,2).dem01);
age = data(:,3).dem02;
lmenis = table2array(data(:,4));
bmi = data(:,5).bmi2;
%
% Read Results Spreadsheet
%
xlsnam = fullfile('Results','mri_fitps','mri_fitps.xlsx');
res = readtable(xlsnam);
%
% Get Variable Names (Column Headers)
%
hdrs = res.Properties.VariableDescriptions;
hdrs = [hdrs(1) {'Sex' 'Age' 'BMI' 'LatMeniscus'} hdrs(2:end)];
%
% Expand Subject Demographics to All Trials
%
subjs = res(:,1).Subject;
%
idx = double(subjs==subjn');
sex = idx*sex;
age = idx*age;
bmi = idx*bmi;
lmenis = idx*lmenis;
%
% Expand Table with Demographics
%
res = addvars(res,sex,age,bmi,lmenis,'After','Subject');
res.Properties.VariableNames = hdrs;
%
% Write Table to MS-Excel Spreadsheet
%
xlsnams = fullfile('Results','mri_fitps','mri_fitps_stat.xlsx');
writetable(res,xlsnams,'WriteMode','replacefile');
%
return