%
% Get Subject Directories
%
rdirs = { '005'   % 1
          '006'   % 2
          '007'   % 3 
          '008'   % 4
%%           '010' % 5     % 010_R_ALL_RHO_Y1_MK.csv: Lateral compartment segmentations only - no medial compartment segmentations.
          '011'   % 6
          '014'   % 7
          '015' % 8     % 015_L_ALL_RHO_Y1_MK.csv: Divot in femoral cartilage. Needs a point in the bone ROI at cartilage divot?
          '017'   % 9
          '018' % 10    % 018_R_ALL_T2S_Y1_EF.csv: Need additional bone points where the cartilage is rough on both the femur and tibia on slice 39 (and others?).
          '020'   % 11
          '021'   % 12
          '023' % 13    % 023_L_ALL_T2S_Y1_AB.csv: Medial condyle bone ROI named "c" instead of "MEDCONDYLE."
          '024'   % 14
          '025'   % 15
          '027'   % 16
          '028' % 17    % 028_L_ALL_RHO_Y1_MK.csv: Trochlea slice 32 needs a point in the bone ROI at cartilage divot.  See Subj028_T1rho_plt.pdf.
          '031'   % 18
          '033' % 19    % 033_R_ALL_RHO_Y1_KP.csv: Tripled tibia cartilage ROIs (LACS) on slices 58, 59, and 60.
          '034' % 20    % 034_L_AX_TIB_T2S_Y1_AB.csv: Too many POSTAXIS slices and no DIST_TIBIA ROI.  Labeling problem?
          '035' % 21    % No axial files!!! (RHO\035_L_AX_FEM*.csv, RHO\035_L_AX_TIB*.csv, RHO\035_R_AX_FEM*.csv, and RHO\035_R_AX_TIB*.csv)
          '037' % 22    % 037_L_ALL_T2S_Y1_EF.csv: On slice 41, the ROI for "MEDCONDYLECART" has a space at the beginning of the ROI name: " MEDCONDYLECART".
          '038'   % 23
          '040' % 24    % 040_L_ALL_T2S_Y1_AB.csv: On slices 6-16 there are two MEDCONDYLE ROIs and no MEDCONDYLECART ROIs.  One needs renaming.
          '042' % 25    % 042_L_ALL_T2S_Y1_EF.csv: On slice 15 (others?), need additional bone points where the cartilage is rough (a couple of divots) on the posterior femur.
          '049'   % 26
          '051' % 27    % 051_R_ALL_T2S_Y1_EF.csv:  Almost all of the femoral cartilage ROIs (LATCONDYLECART, MEDCONDYLECART, and TROCHLEACART) are duplicated.  Please delete the duplicates.
          '053' % 28    % 053_L_ALL_RHO_Y1_MK.csv: Hole in posterior femoral cartilage.  Cartilage segmentation truncated.
          '055'   % 29
%%           '056' % 30    % 056_L_ALL_RHO_Y1_MK.csv: No lateral compartment segmentations.
          '057' % 31    % 057_L_AX_TIB_T2S_Y1_AB.csv: No distal tibia ROI.
          '058' % 32    % 058_R_ALL_T2S_Y1_AB.csv: No LATCONDYLE. Mislabeled?
          '060'   % 33
          '061' % 34    % 061_L_AX_TIB_RHO_Y1_MK.csv: Too many POSTAXIS slices.
%%           '062' % 35    % No left leg RHO files (062_L_ALL_RHO*.csv and 062_L_AX_*.csv files).  Too dark to segment!
          '063' % 36    % 063_L_ALL_T2S_Y1_AB.csv: On slice 38, the posterior tibia bone line reverses on itself. See PTOA_063_L_Y1_T2S_Slice38_TibiaBone.pdf.
          '064' % 37    % 064_L_ALL_RHO_Y1_MK.csv: No LATCONDYLECART. Mislabeled as MEDCONDYLECART?
          '067'   % 38
          '070'   % 39
          '077' % 40    % 077_L_ALL_RHO_Y1_MK.csv: Slice 8 has three "LACS" ROIs which are all exactly the same and all overlapping. 
          '080' % 41    % 080_L_ALL_RHO_Y1_MK.csv: Divot in femoral cartilage. Needs one or two point in the bone ROI at cartilage divot?
          '082'   % 42
          '091' % 43    % 091_L_ALL_T2S_Y1_MK.csv: No femoral segmentations, just tibia?
          '092'   % 44
          '096' };% 45
%
rdirs = rdirs([2:2:6 11]);             % Redo 6, 8, 14, and 21
%
nsubj = size(rdirs,1);
bdir = 'YEAR1';
%
tim = zeros(nsubj,1);   % Time to create masks
%
for kk = 1:nsubj
%
   cd(fullfile(rdirs{kk},bdir));
%
   pwd
%
   tstart = tic;
   seg_prois_b;
   close all;
   clear brois d* f* i* l* m* nm* nr* nsl* p* rdir rm* rs* rt* rx s*;
   clear t t2 t3* tp* v x*;
   tim(kk) = toc(tstart);
%
   cd ..\..
%
end
%
tim
%
return
