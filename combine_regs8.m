%#######################################################################
%
%               * Combine REGistrations Data 8 Program *
%
%          M-File which reads the registration data from different
%     subject MS-Excel spreadsheet files and combines them as different
%     sheets in the output MS-Excel spreadsheet file,
%     Registration_data*.xlsx.
%
%     NOTES:  1.  The program should be run from the top level
%             directory for the PTOA study.
%
%     20-Feb-2025 * Mack Gardner-Morse
%

%#######################################################################
%
% Output MS-Excel Spreadsheet File
%
ddir = fullfile('Results','Registration');       % Data directory
xlsf = 'Registration_data_19Feb2025.xlsx';       % Output spreadsheet file name
xlsf = fullfile(ddir,xlsf);
%
% Subject and Visit Directories
%
sdirs = {'005\YEAR2'
         '006\YEAR2'
         '007\YEAR2'
         '008\YEAR2'
         '015\YEAR1'
         '015\YEAR2'
         '017\YEAR2'
         '020\YEAR2'
         '021\YEAR2'
         '023\YEAR1'
         '023\YEAR2'
         '027\YEAR1'
         '027\YEAR2' };

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