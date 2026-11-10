%#######################################################################
%
%                     * Combine SHeeTS 5 Program *
%
%          M-File which reads the subject data from different sheets in
%     the MS-Excel spreadsheet file, Registration_data_05Aug2024.xlsx,
%     and combines the subject data into one sheet in the output
%     MS-Excel spreadsheet file, Combined_Registration_05Aug2024.xlsx.
%
%     NOTES:  None.
%
%     26-Jul-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Input and Output MS-Excel Spreadsheet Files
%
ddir = fullfile('Results','Registration');       % Data directory
xlsf = 'Registration_data_05Aug2024.xlsx';       % Output spreadsheet file name
xlsf = fullfile(ddir,xlsf);
%
xlsf2 = 'Combined_Registration_05Aug2024.xlsx';  % Output spreadsheet file name
xlsf2 = fullfile(ddir,xlsf2);
%
% Get Sheets in Input Spreadsheet File
%
[~,shts] = xlsfinfo(xlsf);
nsht = length(shts);
%
% Loop through Sheets
%
for k = 1:nsht
%
   sht = shts{k};       % Current sheet
   subj_visit = regexp(sht,'_','split');    % Get subject ID and visit ID
   subj = subj_visit{1};
   visit = subj_visit{2};
%
% Read Subject Data from Input Spreadsheet File
%
   raw = readcell(xlsf,'Sheet',shts{k});    % Read subject data
%
   nr = size(raw,1);
   raw = [cellstr(repmat(subj,nr,1)) cellstr(repmat(visit,nr,1)) raw];
   if k>1
     raw = raw(2:end,:);
   end
%
% Write Subject Data to New Output Spreadsheet File
%
   writecell(raw,xlsf2,'Sheet','Combined','WriteMode','append');
%
end