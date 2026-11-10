%#######################################################################
%
%               * Combine REGistrations Data 4 Program *
%
%          M-File which reads the registration data from different
%     subject MS-Excel spreadsheet files and combines them as different
%     sheets in the output MS-Excel spreadsheet file,
%     Registration_data?.xlsx.
%
%     NOTES:  1.  The program should be run from the top level
%             directory for the PTOA study.
%
%     05-Aug-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Output MS-Excel Spreadsheet File
%
ddir = fullfile('Results','Registration');       % Data directory
xlsf = 'Registration_data_05Aug2024.xlsx';       % Output spreadsheet file name
xlsf = fullfile(ddir,xlsf);
%
% Subject and Visit Directories
%
sdirs = {'015\BASELINE'
         '022\BASELINE'
         '033\BASELINE'
         '037\BASELINE'
         '040\BASELINE'
         '043\BASELINE'
         '047\BASELINE'
         '049\BASELINE'
         '056\BASELINE'
         '058\BASELINE'
         '061\BASELINE' };
nreg = size(sdirs,1);   % Number of registrations
%
% Loop Through Subject and Visit Spreadsheets
%
for kk = 1:nreg
%
   xlsdir = sdirs{kk};
   dirstrs = regexp(xlsdir,'\','split');
   dirstr = dirstrs{end};
%
   xlsnam = fullfile(xlsdir,[dirstr '.xlsx']);
%
   shtnam = replace(xlsdir,filesep,'_');
%
% Read and Write Sheet
%
   raw = readcell(xlsnam);
   writecell(raw,xlsf,'Sheet',shtnam);
%
 end