function bones = rd_prois(fnam)
%RD_PROIS  Reads a PTOA leg sagittal segmentations for the femur and
%          tibia in 2-D pixels and 3-D mm coordinates from a CSV
%          segmentation file.
%
%          For creating regions of interests (ROIs) of knee joint bone
%          and cartilage for the MRI PTOA study.  ROIs are returned in
%          a structure BONES.  BONES(1) is the femur and BONES(2) is
%          the tiba.  If there are bloom artifacts, BONES(3) are
%          regions of bloom artifact.
%
%          BONES.rois(1) is the cartilage, BONES.rios(2) is the bone,
%          and BONES.rois(3) is the bloom if there are bloom artifact
%          segmentations.
%
%          For the tibia and femur, BONES.rois.roi(1) is the lateral
%          compartment and BONES.rois.roi(2) is the medial compartment.
%          For the femur, BONES(1).rois.roi(3) is the femoral trochlea
%          if there are trochlea segmentations.
%
%          BONES = RD_PROIS(FNAM) given the CSV segmentation file name,
%          PNAM, return structure BONES with the femur, tibia, and
%          bloom artifact segmented regions of interest (ROIs).
%
%          NOTES:  1.  M-file rd_roi6hp.m must be in the current
%                  directory or path.
%
%          15-Sep-2022 * Mack Gardner-Morse 
%
%          02-Jul-2024 * Mack Gardner-Morse * Updated function rd_roi6.m
%                                             with rd_roi6hp.m.
%
%          08-Aug-2025 * Mack Gardner-Morse * Traps for segmentations
%                                             without trochlea.
%
%          26-May-2026 * Mack Gardner-Morse * Added bloom artifact ROI
%                                             segmentations.  See line
%                                             63 and lines 86-103.
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
dat = rd_roi6hp(fnam,true);
dat3 = rd_roi6hp(fnam,false);
%
nams = {dat.name}';
data = {dat.data}';     % 2-D pixels
imageno = {dat.imageno}';
%
data3 = {dat3.data}';   % 3-D coordinates in mm
%
% Get Indices to Bloom Artifact, Femur, and Tibia Data
%
idb = contains(nams,'BLOOM','IgnoreCase',true);  % Index to bloom   
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
idt = ~idf&~idb;        % Index to tibial segmentations
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
% Put Bloom Artifact ROIs into Structures
%
if any(idb)
%
  roi(1).data = data{idb}';            % Bloom
  roi(1).data3 = data3{idb}';
  roi(1).imageno = imageno{idb};
  roi(1).name = nams{idb};
%
  rois(1).name = 'Bloom';
  rois(1).roi = roi;
  rois(1).slice = [imageno{idb}]';
%
  bones(3).name = 'Bloom';
  bones(3).rois = rois;
%
end
%
% Put Tibia ROIs into Structures
%
idtmb = idtm&idtb;      % Tibial medial compartment bone
roi(2).data = data{idtmb}';            % Medial compartment bone
roi(2).data3 = data3{idtmb}';
roi(2).imageno = imageno{idtmb};
roi(2).name = nams{idtmb};
%
idtlb = idtl&idtb;      % Tibial lateral compartment bone
roi(1).data = data{idtlb}';            % Lateral compartment bone
roi(1).data3 = data3{idtlb}';
roi(1).imageno = imageno{idtlb};
roi(1).name = nams{idtlb};
%
rois(2).name = 'Bone';
rois(2).roi = roi;
rois(2).slice = [imageno{(idtl|idtm)&idtb}]';
%
idtmc = idtm&idtc;      % Tibial medial compartment cartilage
roi(2).data = data{idtmc}';            % Medial compartment cartilage
roi(2).data3 = data3{idtmc}';
roi(2).imageno = imageno{idtmc};
roi(2).name = nams{idtmc};
%
idtlc = idtl&idtc;      % Tibial lateral compartment cartilage
roi(1).data = data{idtlc}';            % Lateral compartment cartilage
roi(1).data3 = data3{idtlc}';
roi(1).imageno = imageno{idtlc};
roi(1).name = nams{idtlc};
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
if any(idft)
  idftb = idft&idfb;    % Femoral trochlea bone
  roi(3).data = data{idftb}';          % Trochlea bone
  roi(3).data3 = data3{idftb}';
  roi(3).imageno = imageno{idftb};
  roi(3).name = nams{idft&idfb};
end
%
idfmb = idfm&idfb;      % Femoral medial compartment bone
roi(2).data = data{idfmb}';            % Medial compartment bone
roi(2).data3 = data3{idfmb}';
roi(2).imageno = imageno{idfmb};
roi(2).name = nams{idfmb};
%
idflb = idfl&idfb;      % Femoral lateral compartment bone
roi(1).data = data{idflb}';            % Lateral compartment bone
roi(1).data3 = data3{idflb}';
roi(1).imageno = imageno{idflb};
roi(1).name = nams{idflb};
%
rois(2).name = 'Bone';
rois(2).roi = roi;
rois(2).slice = [imageno{(idfl|idfm|idft)&idfb}]';
%
if any(idft)
  idftc = idft&idfc;    % Femoral trochlea cartilage
  roi(3).data = data{idftc}';          % Trochlea cartilage
  roi(3).data3 = data3{idftc}';
  roi(3).imageno = imageno{idftc};
  roi(3).name = nams{idftc};
end
%
idfmc = idfm&idfc;      % Femoral medial compartment cartilage
roi(2).data = data{idfmc}';            % Medial compartment cartilage
roi(2).data3 = data3{idfmc}';
roi(2).imageno = imageno{idfmc};
roi(2).name = nams{idfmc};
%
idflc = idfl&idfc;      % Femoral lateral compartment cartilage
roi(1).data = data{idflc}';            % Lateral compartment cartilage
roi(1).data3 = data3{idflc}';
roi(1).imageno = imageno{idflc};
roi(1).name = nams{idflc};
%
rois(1).name = 'Cartilage';
rois(1).roi = roi;
rois(1).slice = [imageno{(idfl|idfm|idft)&idfc}]';
%
bones(1).name = 'Femur';
bones(1).rois = rois;
%
return