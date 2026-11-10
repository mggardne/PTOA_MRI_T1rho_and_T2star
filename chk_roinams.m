%
% ROI Names
%
rnams = {'LACS'
         'LATCONDYLE'
         'LATCONDYLECART'
         'LTCS'
         'MACS'
         'MEDCONDYLE'
         'MEDCONDYLECART'
         'MTCS'
         'TROCHLEA'
         'TROCHLEACART' };
%
% T1rho and T2* Directories
%
rdirs = {'RHO'          % Directory for T1rho segmentations
         'T2S' };       % Directory for T2* segmentations
%
% Subject and Visit Directory with CSV Files as a String
%
dirstr = split(pwd,filesep);
%
% Get Subject
%
subj = dirstr{end-1};
%
% Loop through ROI Directories
%
for k = 1:2
%
   rdir = rdirs{k};     % Directory for segmentations
%
% Get CSV Files
%
   d = dir(fullfile(rdir,[subj '*_ALL_*.csv']));
   cnams = {d.name}';      % CSV file names
%
   ncsv = size(cnams,1);
%
% Loop through CSV Files
%
   for l = 1:ncsv
      dat = rd_roi6hp(fullfile(rdir,cnams{l}));
      nnam = size(dat,1);
      rois = {dat.name}';
%
% Check for Bloom ROI
%
      ibloom = any(startsWith(rois,'BLOOM','IgnoreCase',true));
      if ibloom
        rnams_chk = {'BLOOM' rnams{:}}';
        ncorrect = 11;
      else
        rnams_chk = rnams;
        ncorrect = 10;
      end
%
% Check Number of ROIs
%
      if nnam~=ncorrect
        fprintf(1,'\n%s has %i ROI names!\n',cnams{l},nnam);
      end
%
% Check for Missing ROIs
%
      miss = setdiff(rnams_chk,rois);
      if ~isempty(miss)
        fprintf(1,'\n%s is missing the following ROI name(s):\n', ...
                cnams{l});
        for m = 1:length(miss)
           fprintf(1,'%s\n',miss{m});
        end
        fprintf(1,'\n');
      end
%
% Check for Extra ROIs
%
      extra = setdiff(rois,rnams_chk);
      if ~isempty(extra)
        fprintf(1,'\n%s has the following extra ROI name(s):\n', ...
                cnams{l});
        for m = 1:length(extra)
           fprintf(1,'%s\n',extra{m});
        end
        fprintf(1,'\n');
      end
   end
%
end
%
return