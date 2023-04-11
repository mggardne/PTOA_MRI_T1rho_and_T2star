%#######################################################################
%
%               * Combine REGistrations Data 2 Program *
%
%          M-File which reads the registration data from different
%     subject MS-Excel spreadsheet files and combines them as different
%     sheets in the output MS-Excel spreadsheet file,
%     Registration_data?.xlsx.
%
%     NOTES:  1.  The program should be run from the top level
%             directory for the PTOA study:
%             \MRI Postprocessing\
%
%     08-Sep-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Output MS-Excel Spreadsheet File
%
ddir = fullfile('Results','Posters');  % Data directory
xlsf = 'Registration_data.xlsx';       % Output spreadsheet file name
xlsf = fullfile(ddir,xlsf);
%
% Get Subject Directories and Visit Subdirectories
%
subjn = [3;7;12;17;18;19;20;21;24;30;31;32;34;35;38;46];
nsubj = size(subjn,1);
%
sdirs = int2str(subjn);
sdirs = [repmat('0',nsubj,1) sdirs];
sdirs = strrep(cellstr(sdirs),' ','0');
%
% Loop Through Subject and Visit Spreadsheets
%
for kk = 1:nsubj
%
   xlsdir = sdirs{kk};
   dirstr = xlsdir;
%
   xlsnam = fullfile(xlsdir,[dirstr '.xlsx']);
%
   shtnam = dirstr;
%
% Read and Write Sheet
%
   raw = readcell(xlsnam);
   writecell(raw,xlsf,'Sheet',shtnam);
%
 end