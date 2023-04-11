%#######################################################################
%
%                     * PLoT CSV Files 6 Program *
%
%          M-File which reads the cartilage and bone digitization CSV
%     files and plots the data for visual verification.  The program
%     assumes the femur has three sets of digitizations (lateral,
%     medial and trochlea), the tibia has two sets of digitizations
%     (lateral and medial) and the patella has one set of digitizations.
%     Due to the angle of the knee alignment with the scanner, the
%     patella data is ignored.
%
%     NOTES:  1.  Matlab M-files plt_datsl.m and rd_roi6.m must be in
%             the current directory or path.
%
%             2.  Plots are output to PS files:  dig_plt6*.ps
%
%             3.  For the MRI reliability study subjects with left and
%             right knees in the directory.
%             The CSV file names must contain the following key:
%             Side:  _L_ for left knees or _R_ for right knees
%
%     06-Sep-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Directories with Digitization CSV Files
%
dirstr = split(pwd,filesep);
dirstr = dirstr{end};                  % Subject identifier
%
mdirs = ['RHO'; 'T2S'];                % MRI type directories
mtxts = {'T1\rho'; 'T2*'};
nmdirs = size(mdirs,1);
%
% Side Code in File Names
%
sides = ['_L_'; '_R_'];                % L - left and R - right
stxts = {'Left'; 'Right'};
%
% Colors for Different Segmentations
%
% map = colormap;         % Default color map
map = jet;              % Jet color map
midx = 1:26:256;        % Index into color map
midx = midx([2 1 3:7 10 8 9]);
%
% Loop Through Directories and Leg, and Plot Digitizations
%
psfile = ['dig_plt6_' dirstr '.ps'];   % PS file name
%
for l = 1:nmdirs
%
   mdir = mdirs(l,:);
%
   for m = 1:2       % Side (L/R)
%
      side = sides(m,:);
%
% Get CSV File
%
      csv_file = [dirstr side 'ALL_' mdir '_*.csv'];
      dnam = dir(fullfile(mdir,csv_file));
      dnam = {dnam.name}';
%
% Get the Digitization Data
%
      dat = rd_roi6(fullfile(mdir,dnam{1}));
%
      nams = {dat.name}';
      data = {dat.data}';
%
      nd = size(nams,1);
%
% Set Up Figure
%
      figure;
      orient landscape;
%
      h = cell(nd,1);
      hl = gobjects(nd,1);
%
% Plot All the Segmentations
%
      for k = 1:nd
         datc = data{k};
         h{k} = plt_datsl(datc,'k.-',0.5,map(midx(k),:));
         hl(k) = h{k}(1);
      end
%
      legend(hl,nams);
      axis equal;
%
      st = {['Subject ' dirstr]; mtxts{l}; [stxts{m} ' Leg']};
      title(st,'FontSize',16,'FontWeight','bold');
      orient landscape;
      if m==1&&l==1
        print('-dpsc2','-r600','-fillpage',psfile); 
      else
        print('-dpsc2','-r600','-fillpage','-append',psfile); 
      end
   end
end
%
pause;
%
close all;
%
return