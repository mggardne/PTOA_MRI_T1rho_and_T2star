function bones = rd_prois(fnam)
%RD_PROIS  Reads a PTOA leg sagittal segmentations for the femur and
%          tibia in 2-D pixels and 3-D mm coordinates from a CSV
%          segmentation file.
%
%          For creating regions of interests (ROIs) of knee joint bone
%          and cartilage for the MRI PTOA study.  ROIs are returned in
%          a structure BONES.  BONES(1) is the femur and BONES(2) is
%          the tiba.  BONES.rois(1) is the cartilage and BONES.rios(2)
%          is the bone.  BONES.rois.roi(1) is the lateral compartment and
%          and BONES.rois.roi(2) is the medial compartment.  For the
%          femur, BONES(1).rois.roi(3) is the femoral trochlea.
%
%          BONES = RD_PROIS(FNAM) given the CSV segmentation file name,
%          PNAM, return structure BONES with the femur and tibia
%          segmented regions of interest (ROIs).
%
%          NOTES:  1.  M-file rd_roi6.m must be in the current directory
%                  or path.
%
%          15-Sep-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check for Input CSV Segmentation File Name
%
if (nargin<1)
  error([' *** ERROR in RD_PROIS:  An input CSV segmentation file', ...
         ' name is required!']);
end
%
% Read Segmentations
%
dat = rd_roi6(fnam,true);
dat3 = rd_roi6(fnam,false);
%
nams = {dat.name}';
data = {dat.data}';     % 2-D pixels
imageno = {dat.imageno}';
%
data3 = {dat3.data}';   % 3-D coordinates in mm
%
% Get Indices to Femur and Tibia Data
%
idc = contains(nams,'COND','IgnoreCase',true);        % Index to condyles
idft = contains(nams,'TROCH','IgnoreCase',true);      % Index to trochlea
idf = idft|idc;         % Index to femoral segmentations
idcl = startsWith(nams(idc),'LAT','IgnoreCase',true); % Lateral condyles
idfl = idc;
idfl(idc) = idcl;       % Index to lateral condyles
idfm = idc;
idfm(idc) = ~idcl;      % Index to medial condyles
idfc = endsWith(nams,'CART','IgnoreCase',true);       % Index to femoral cartilage
idfb = idf&~idfc;       % Index to femoral bone
%
idt = ~idf;             % Index to tibial segmentations
ida = contains(nams(idt),'A','IgnoreCase',true);      % Index to cartilage
idtc = idt;
idtc(idt) = ida;        % Index to tibial cartilage
idtb = idt;
idtb(idt) = ~ida;       % Index to tibial bone
idl = startsWith(nams(idt),'L','IgnoreCase',true);
idtl = idt;
idtl(idt) = idl;        % Index to lateral tibia
idtm = idt&~idtl;       % Index to medial tibia
%
% Put Tibia ROIs into Structures
%
roi(2).data = data{idtm&idtb}';        % Medial compartment bone
roi(2).data3 = data3{idtm&idtb}';
roi(2).imageno = imageno{idtm&idtb};
roi(2).name = nams{idtm&idtb};
%
roi(1).data = data{idtl&idtb}';        % Lateral compartment bone
roi(1).data3 = data3{idtl&idtb}';
roi(1).imageno = imageno{idtl&idtb};
roi(1).name = nams{idtl&idtb};
%
rois(2).name = 'Bone';
rois(2).roi = roi;
rois(2).slice = [imageno{(idtl|idtm)&idtb}]';
%
roi(2).data = data{idtm&idtc}';        % Medial compartment cartilage
roi(2).data3 = data3{idtm&idtc}';
roi(2).imageno = imageno{idtm&idtc};
roi(2).name = nams{idtm&idtc};
%
roi(1).data = data{idtl&idtc}';        % Lateral compartment cartilage
roi(1).data3 = data3{idtl&idtc}';
roi(1).imageno = imageno{idtl&idtc};
roi(1).name = nams{idtl&idtc};
%
rois(1).name = 'Cartilage';
rois(1).roi = roi;
rois(1).slice = [imageno{(idtl|idtm)&idtc}]';
%
bones(2).name = 'Tibia';
bones(2).rois = rois;
%
% Put Femur ROIs into Structures
%
roi(3).data = data{idft&idfb}';        % Trochlea bone
roi(3).data3 = data3{idft&idfb}';
roi(3).imageno = imageno{idft&idfb};
roi(3).name = nams{idft&idfb};
%
roi(2).data = data{idfm&idfb}';        % Medial compartment bone
roi(2).data3 = data3{idfm&idfb}';
roi(2).imageno = imageno{idfm&idfb};
roi(2).name = nams{idfm&idfb};
%
roi(1).data = data{idfl&idfb}';        % Lateral compartment bone
roi(1).data3 = data3{idfl&idfb}';
roi(1).imageno = imageno{idfl&idfb};
roi(1).name = nams{idfl&idfb};
%
rois(2).name = 'Bone';
rois(2).roi = roi;
rois(2).slice = [imageno{(idfl|idfm|idft)&idfb}]';
%
roi(3).data = data{idft&idfc}';        % Trochlea cartilage
roi(3).data3 = data3{idft&idfc}';
roi(3).imageno = imageno{idft&idfc};
roi(3).name = nams{idft&idfc};
%
roi(2).data = data{idfm&idfc}';        % Medial compartment cartilage
roi(2).data3 = data3{idfm&idfc}';
roi(2).imageno = imageno{idfm&idfc};
roi(2).name = nams{idfm&idfc};
%
roi(1).data = data{idfl&idfc}';        % Lateral compartment cartilage
roi(1).data3 = data3{idfl&idfc}';
roi(1).imageno = imageno{idfl&idfc};
roi(1).name = nams{idfl&idfc};
%
rois(1).name = 'Cartilage';
rois(1).roi = roi;
rois(1).slice = [imageno{(idfl|idfm|idft)&idfc}]';
%
bones(1).name = 'Femur';
bones(1).rois = rois;
%
return