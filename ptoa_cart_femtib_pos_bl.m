%#######################################################################
%
%       * PTOA CARTilage FEMurTIBia POSition BaseLine Program *
%
%          M-File which reads the PTOA femur and tibia relative position
%     data from the MS-Excel spreadsheet file in the Results\FemTib_Pos
%     subdirectory:  cthick_tibfem_pos4_BL.xlsx.
%
%          The PTOA cartilage thickness results are read and the
%     relative positions are incorporated into the results.  The
%     expanded results are written to the MS-Excel spreadsheets:
%     Femur Thickness Averages baseline_rp.xlsx and
%     Tibia Thickness Averages baseline_rp.xlsx in the
%     Results\Osee_CartThick\ subdirectory.
%
%     NOTES:  1.  Assumes the default directory is PTOA\Results.
%
%     13-Nov-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Read Relative Position Data for Baseline PTOA Subjects from MS-Excel
% Spreadsheet File
%
xlsnm1 = 'cthick_tibfem_pos4_BL.xlsx'; % Tibia wrt the femur
% xlsnm2 = 'cthick_femtib_pos4_BL.xlsx';
pos_dat = readtable(fullfile('FemTib_Pos',xlsnm1));
%
pos_subj = pos_dat.Subject;  % Subject ID as a number
pos_res = pos_dat.Result;    % Analysis result (0 = T1rho and 1 = T2*)
pos_leg = pos_dat.Leg;       % Leg (0 = Left and 1 = Right)
%
idx = ~logical(pos_res);     % Index to T1rho data
%
pos_subj = pos_subj(idx);    % T1rho subjects
pos_leg = pos_leg(idx);      % T1rho legs
%
pos = pos_dat{idx,7:12};     % Relative position data for T1rho
%
% Read Results Spreadsheets
%
xlsnamf = fullfile('Osee_CartThick', ...
                   'Femur Thickness Averages baseline.xlsx');
xlsnamt = fullfile('Osee_CartThick', ...
                   'Tibial Thickness Averages baseline.xlsx');
res_datf = readtable(xlsnamf);
res_datt = readtable(xlsnamt);
%
% Get Variable Names (Column Headers)
%
% hdrs = res_dat.Properties.VariableNames;
% hdrs = [hdrs {'ap' 'ml' 'si' 'ir' 'fe' 'vv'}];
hdrs = {'ap' 'ml' 'si' 'ir' 'fe' 'vv'};
%
% Expand Relative Position Data to All ROIs
%
res_subjf = res_datf.Subject;          % Subject ID as a number
res_legf = res_datf.Leg;               % Leg (0 = Left and 1 = Right)
%
res_subjt = res_datt.Subject;          % Subject ID as a number
res_legt = res_datt.Leg;               % Leg (0 = Left and 1 = Right)
%
id1f = res_subjf==pos_subj';
id2f = res_legf==pos_leg';
idf = id1f&id2f;        % Matrix relating position data to ROIs
%
res_posf = idf*pos;     % Expanded relative position data for the femur
%
id1t = res_subjt==pos_subj';
id2t = res_legt==pos_leg';
idt = id1t&id2t;        % Matrix relating position data to ROIs
%
res_post = idt*pos;     % Expanded relative position data for the tibia
%
% Expand Table with Relative Positions
%
rf = array2table(res_posf);
rf.Properties.VariableNames = hdrs;
tf = [res_datf rf];
%
rt = array2table(res_post);
rt.Properties.VariableNames = hdrs;
tt = [res_datt rt];
%
% Write Table to MS-Excel Spreadsheet
%
idot = strfind(xlsnamf,'.');
idot = idot(end);
xlsnamof = [xlsnamf(1:idot-1) '_rp' xlsnamf(idot:end)];
writetable(tf,xlsnamof,'WriteMode','replacefile');
%
idot = strfind(xlsnamt,'.');
idot = idot(end);
xlsnamot = [xlsnamt(1:idot-1) '_rp' xlsnamt(idot:end)];
writetable(tt,xlsnamot,'WriteMode','replacefile');
%
return