function bones = rd_rois(rdir,leg,ld,itroch,irho)
%RD_ROIS   Reads a particular loaded or unloaded leg sagittal
%          segmentation for the femur and tibia in pixels from CSV
%          files in a directory.
%
%          For creating regions of interest ROIs of knee joint cartilage
%          for the MRI reliability study.  ROIs are returned in a
%          structure BONES.  BONES(1) is the femur and BONES(2) is the
%          tiba.  BONES.rois(1) is the cartilage and BONES.rios(2) is
%          the bone.  BONES.rois.roi(1) is the lateral compartment and
%          BONES.rois.roi(2) is the medial compartment.
%
%          BONES = RD_ROIS(RDIR,LEG,LD) given the directory name in the
%          string, RDIR, either the character 'L' or 'R' for the left
%          or right leg in LEG, and either 'LD' or 'UL' for loaded or
%          unloaded condition in LD, return structure BONES with the
%          femur and tibia segmented regions of interest (ROIs).
%
%          BONES = RD_ROIS(RDIR,LEG,LD,ITROCH) given the
%          logical ITROCH, trochlear ROIs are returned in the femur
%          structure BONES(1).
%
%          BONES = RD_ROIS(RDIR,LEG,LD,ITROCH,IRHO) given the
%          integer IRHO, checks for 'imageno' greater than 96 and
%          subtracts one (1) from 'imageno', divides by IRHO, and adds
%          one (1).  This is to account for the digitization on the
%          first spin lock time.
%
%          NOTES:  1.  M-file rd_roi6.m must be in the current directory
%                  or path.
%
%                  2.  Femur and tibia segmentations must be in sub-
%                  directories 'Femur' and 'Tibia', respectively under
%                  segmentation directory RDIR.
%
%                  3.  For sagittal segmentations only.  CSV file names
%                  must contain the capital letters 'SAG'.
%
%                  4.  Trochlear ROIs are removed by default.
%
%                  5.  CSV files with case insensitive "dup" in the
%                  file names are ignored as duplicate files.
%
%                  6.  Cartilage (SAGAR) CSV files with "_RO" in the
%                  file names are used in place of the same files
%                  without "_RO".
%
%                  7.  Trochlea option is NOT tested.
%
%          20-Jan-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in RD_ROIS:  Three inputs are required!');
end
%
if (nargin<4)
  itroch = false;
end
%
if (nargin<5)
  irho = 1;
end
%
% Loop through the Femur and Tibia Bones
%
bdirs = ['Femur'; 'Tibia'];
nb = size(bdirs,1);
%
bones = struct;
%
for l = 1:nb
%
   rnams = ['*' leg '*SAG*' ld '*.csv'];
   rnams = dir(fullfile(rdir,bdirs(l,:),rnams));
   rnams = {rnams.name}';
   idx = ~contains(rnams,'MGG','IgnoreCase',true);
   rnams = rnams(idx);
   idx = ~contains(rnams,'dup','IgnoreCase',true);
   rnams = sort(rnams(idx));           % Cartilage first, bone last
   idx = contains(rnams,'_RO');        % Check for _RO files
   if any(idx)
     idc = contains(rnams,'SAGAR');
     idx = idx|~idc;
     rnams = rnams(idx);
   end
   nrfiles = size(rnams,1);
   if nrfiles>2
     fprintf(1,'  %s\n',rnams{:});
     error(' *** ERROR in RD_ROIS:  Too many CSV files found!');
   end
%
% Loop through ROI Files
% Note trochlear slices are removed by default.
%
   rois = struct;
%
   for k = 1:nrfiles
%
      rnam = fullfile(rdir,bdirs(l,:),rnams{k});
%
      rois(k).name = rnams{k};
      rois(k).roi = rd_roi6(rnam,true);     % Read pixel coordinates
      if l==1&&~itroch
        rois(k).roi = rois(k).roi(1:2);     % Remove trochlear slices
      end
      rois(k).slice = [rois(k).roi.imageno]';
      if any(rois(k).slice>96)&&irho>1
        rois(k).slice = (rois(k).slice-1)./irho+1;
        rois(k).roi(1).imageno = (rois(k).roi(1).imageno-1)./irho+1;
        rois(k).roi(2).imageno = (rois(k).roi(2).imageno-1)./irho+1;
        if l==1&&itroch
          rois(k).roi(3).imageno = (rois(k).roi(3).imageno-1)./irho+1;
        end
      end
%
   end
%
   bones(l).name = bdirs(l,:);
   bones(l).rois = rois;
%
end
%
return