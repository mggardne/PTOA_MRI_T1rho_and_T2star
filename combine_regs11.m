%#######################################################################
%
%              * Combine REGistrations Data 11 Program *
%
%          M-File which reads the registration data from different
%     subject MS-Excel spreadsheet files and combines them as different
%     sheets in the output MS-Excel spreadsheet file,
%     Registration_data_*.xlsx, where * is the current date.
%
%     NOTES:  1.  The program should be run from the top level
%             directory for the PTOA study.
%
%     25-Jul-2025 * Mack Gardner-Morse
%

%#######################################################################
%
% Output MS-Excel Spreadsheet File
%
dat = char(datetime('today'),'ddMMMyyyy');       % Date
ddir = fullfile('Results','Registration');       % Data directory
xlsf = ['Registration_data_' dat '.xlsx'];       % Output spreadsheet file name
xlsf = fullfile(ddir,xlsf);
%
% Subject and Visit Directories
%
sdirs = {'011\YEAR1'
         '049\YEAR1'
         '056\YEAR1'
         '058\YEAR1' };

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