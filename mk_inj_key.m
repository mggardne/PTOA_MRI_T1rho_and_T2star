%#######################################################################
%
%                     * MaKe INJury KEY Program *
%
%          M-File which reads data from a REDCap injured leg CSV data
%     file and a PTOA MS-Excel results spreadsheet.  The injury data is
%     used to create an injury key column for the PTOA MS-Excel results
%     spreadsheet. 
%
%     NOTES:  1.  The  REDCap injured leg CSV file:
%         PTOAStudyTheAdaptive-InjuredLimb4723_DATA_2023-04-07_1225.csv
%             must be in the current directory or path.
%
%             2.  The PTOA MS-Excel results spreadsheet, mri_fits.xlsx
%             must be in the directory Results\mri_fitps\.
%
%     07-Apr-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% REDCap Injured Leg Data File
%
csvnm = 'PTOAStudyTheAdaptive-InjuredLimb4723_DATA_2023-04-07_1225.csv';
%
% Get Subject Numbers
%
snams = readmatrix(csvnm,'Range','A2:A21','OutputType','char');
snams = char(snams);
snums = str2num(snams(:,5:7));
%
% Get Injury Key
% Injured Right = 1, Injured Left = 2
%
injkey = readmatrix(csvnm,'Range','B2:B21');
%
% Get PTOA Results Subject and Leg Data
%
resdir = 'Results\mri_fitps\';
%
resfile = fullfile(resdir,'mri_fitps.xlsx');
%
subjleg = readmatrix(resfile,'Range','A2');
subjleg = subjleg(:,1:3);
nrows = size(subjleg,1);
%
% Loop through Subjects to Generate Injury Key
%
inj = zeros(nrows,1);
%
for k = 1:20
%
   ks = snums(k);
   ids = subjleg(:,1)==ks;
%
   if injkey(k)==1
     idl = subjleg(ids,3)==1;
   else
     idl = subjleg(ids,3)==0;
   end
%
   inj(ids) = double(idl);
%
end
%
% Write Injury Key Column to Spreadsheet
%
writematrix(inj,resfile,'Range','D2','AutoFitWidth',false);
%
return