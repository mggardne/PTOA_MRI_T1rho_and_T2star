%#######################################################################
%
%             * PTOA FEMurTIBia POSition Year 1 Program *
%
%          M-File which reads the PTOA femur and tibia relative position
%     data from the MS-Excel spreadsheet file in the Results\FemTib_Pos
%     subdirectory:  cthick_tibfem_pos4_Y1.xlsx.
%
%          The PTOA T1rho and T2* results are read and the relative
%     positions are incorporated into the results.  The expanded results
%     are written to the MS-Excel spreadsheet:  mri_fitpsvy1r2rp.xlsx
%     in the Results\mri_fitpsvy1r2\ subdirectory.
%
%     NOTES:  1.  Assumes the default directory is PTOA\Results.
%
%     25-Nov-2025 * Mack Gardner-Morse
%

%#######################################################################
%
% Read Relative Position Data for Year 1 PTOA Subjects from MS-Excel
% Spreadsheet File
%
xlsnm1 = 'cthick_tibfem_pos4_Y1.xlsx'; % Tibia wrt the femur
% xlsnm2 = 'cthick_femtib_pos4_Y1.xlsx';
pos_dat = readtable(fullfile('FemTib_Pos',xlsnm1));
%
pos_subj = pos_dat.Subject;  % Subject ID as a number
pos_res = pos_dat.Result;    % Analysis result (0 = T1rho and 1 = T2*)
pos_leg = pos_dat.Leg;       % Leg (0 = Left and 1 = Right)
%
pos = pos_dat{:,7:12};  % Relative position data
%
% Read Results Spreadsheet
%
xlsnam = fullfile('mri_fitpsvy1r2','mri_fitpsvy1r2.xlsx');
res_dat = readtable(xlsnam);
%
% Get Variable Names (Column Headers)
%
% hdrs = res_dat.Properties.VariableNames;
% hdrs = [hdrs {'ap' 'ml' 'si' 'ir' 'fe' 'vv'}];
hdrs = {'ap' 'ml' 'si' 'ir' 'fe' 'vv'};
%
% Expand Relative Position Data to All ROIs
%
res_subj = res_dat.Subject;  % Subject ID as a number
res_res = res_dat.Result;    % Analysis result (0 = T1rho and 1 = T2*)
res_leg = res_dat.Leg;       % Leg (0 = Left and 1 = Right)
%
id1 = res_subj==pos_subj';
id2 = res_res==pos_res';
id3 = res_leg==pos_leg';
id = id1&id2&id3;       % Matrix relating position data to ROIs
%
res_pos = id*pos;       % Expanded relative position data
%
% Expand Table with Relative Positions
%
r = array2table(res_pos);
r.Properties.VariableNames = hdrs;
t = [res_dat r];
%
% Write Table to MS-Excel Spreadsheet
%
idot = strfind(xlsnam,'.');
idot = idot(end);
xlsnamo = [xlsnam(1:idot-1) 'rp' xlsnam(idot:end)];
writetable(t,xlsnamo,'WriteMode','replacefile');
%
return