function [maskf,maskt,ibone,f,t,f3,t3,maskb,bsl,b,b3] = mk_lay_msks( ...
                             brois,rsl,nrsl,iszs,dist,scal,v,fps,sn,leg)
%MK_LAY_MSKS  Takes the PTOA knee sagittal segmentations for the femur
%          and tibia and creates image masks for deep and superficial
%          cartilage layers.
%
%          [MASKF,MASKT] = MK_LAY_MSKS(BROIS,RSL,NRSL,ISZS) Given a
%          structure, BROIS, with the femur and tibia segmentations,
%          an array of slices with cartilage and bone segmentations,
%          RSL, and the number of slices in RSL, NRSL, returns the
%          femoral mask, MASKF, and tibial mask, MASKT.  The first
%          column in the logical masks is the mask for the image, the
%          second column is the cartilage layers (1 - deep, and 2 -
%          superficial), and the third column is the slices in RSL.
%
%          [MASKF,MASKT] = MK_LAY_MSKS(BROIS,RSL,NRSL,ISZS,DIST,SCAL)
%          Midline points must be within a distance, DIST, of the
%          second line.  The default value is Inf (infinite).  A one or
%          two element scale, SCAL is used to scale the X and Y
%          coordinates before comparison with the distance DIST for
%          midline points.  The scale is usually the pixel size.  See
%          midline.m.
%
%          [MASKF,MASKT,IBONE] = MK_LAY_MSKS(BROIS,RSL,NRSL,ISZS)
%          Returns the logical matrix, IBONE, with the slices in RSL in
%          rows and the first column is the femur and the second column
%          is the tibia.  IBONE is true if the bone is on a slice.
%
%          [MASKF,MASKT,IBONE,F,T] = MK_LAY_MSKS(BROIS,RSL,NRSL,ISZS)
%          Returns the cell arrays, F, and T, with the femoral (F) and
%          tibial (T) segmentations.  The first row of the cell array
%          is the cartilage segmentation, and the second row is the
%          bone segmentations.  The columns are the slices in RSL.
%
%          [MASKF,MASKT,IBONE,F,T,F3,T3] = MK_LAY_MSKS(BROIS,RSL,NRSL,
%          ISZS) Returns the cell arrays, F3, and T3, with the femoral
%          (F3) and tibial (T3) three dimensional (3D) segmentations in
%          the MRI coordinate system.  The first row of the cell array
%          is the cartilage segmentation, and the second row is the
%          bone segmentations.  The columns are the slices in RSL.
%
%          [MASKF,MASKT,IBONE,F,T,F3,T3,MASKB,BSL,B,B3] = MK_LAY_MSKS(
%          BROIS,RSL,NRSL,ISZS) If the structure, BROIS, contains
%          segmentations of bloom artifacts, returns the two-dimensional
%          bloom mask, MASKB. The first column in the logical masks for
%          the image, the second column is the slices in RSL.  The cell
%          arrays, B, and B3, contain the two-dimensional bloom (B) and
%          three-dimensional bloom (B3) segmentations.  The columns are
%          the slices in BSL.
%
%          If the structure, BROIS, does not contain bloom artifact
%          segmentations, the returned bloom mask, MASKB, is a logical
%          true array.  The returned BSL, B, and B3 are empty arrays.
%
%          [MASKF,MASKT] = MK_LAY_MSKS(BROIS,RSL,NRSL,ISZS,DIST,SCAL,
%          V,FPS,SN,LEG) Given the image matrix (third dimension is
%          the slices), V, text for the start of the output file names,
%          FPS, text for the subject number, SN, and a character for
%          the leg (L - left and R - right), LEG, plots the
%          segmentations and masks.  FPS is required to generate
%          Postscript (PS) files. SN and LEG are optional texts that
%          are used to generate plot titles.
%
%          NOTES:  1.  Image matrix V is required to generate plots.
%
%                  2.  FPS is required to generate Postscript (PS)
%                  files.  The plot file names are FPS_rois1.ps and
%                  FPS_rois2.ps.
%
%                  3.  FPS, SN, and LEG are used to title the plots.
%                  Any underscores in FPS are replaced with spaces.
%
%                  4.  The M-files cr_mask2.m, cr_mask2f.m, fix_gap.m,
%                  in_tri2d.m, lsect2.m, lsect2a.m, lsect3.m, lsect4.m,
%                  lsect5.m, midline.m, mk2_tri_2d.m, mk2_tri_2df.m,
%                  and near2.m must be in the current directory or
%                  path.
%
%          09-Dec-2022 * Mack Gardner-Morse
%
%          08-Aug-2025 * Mack Gardner-Morse * Checks for trochlea in
%                                             input BROIS variable.
%
%          04-Jun-2026 * Mack Gardner-Morse * Incorporated bloom
%                                             artifact segmentations to
%                                             create bloom masks.
%

%#######################################################################
%
% Check Number of Inputs
%
if (nargin<6)
  error([' *** ERROR in MK_LAY_MSKS:  At least six (6) inputs are', ...
         ' required!']);
end
%
if (nargin>6)
  iplt = true;
else
  iplt = false;
end
%
if (nargin>7)
  iprt = true;
else
  iprt = false;
end
if isempty(fps)
  fps = 'plt';
end
%
if iprt
  pnam1 = [fps '_ROIs1.ps'];           % ROI lines print file name
  pnam2 = [fps '_ROIs2.ps'];           % ROI areas print file name
%
  ttxt = strrep(fps,'_',' ');
%  
  if exist('sn','var')&&~isempty(sn)
    ttxt = [ttxt ', Subject ' sn];
  end
%  
  if exist('leg','var')
    if strcmpi(leg,'L')
      ttxt = [ttxt ', Left Leg'];
    else
      ttxt = [ttxt ', Right Leg'];
    end
  end
end
%
% Check Inputs
%
if ~isstruct(brois)
  error([' *** ERROR in MK_LAY_MSKS:  First input must be a', ...
         ' structure with the segmentations!']);
end
%
% Setup Color Maps, Line Types and Color, and Legend Strings
%
if iplt
  gmap = gray(128);     % Gray color map for not cartilage
  jmap = jet(128);      % Jet color map for cartilage
  cmap = [gmap; jmap];
%
  lt = ['g.-'; 'b.-'; 'c.-'; 'r.-'; 'y.-'; 'm.-'];    % Line color and type
  legds = ['FEM_CART'; 'FEM_BONE'; 'TIB_CART'; 'TIB_BONE'; 'BLOOM   '];
end
%
% Initialize ROI Coordinate Arrays
%
npxt = prod(iszs);                     % Total number of pixels in a slice
%
if size(brois,2)>2
  ibloom = true;
  bsl = brois(3).rois.slice;           % Slices with bloom (not unique)
  nb = size(bsl,1);     % Number of blooms
  b = cell(1,nb);       % 2D bloom coordinates
  b3 = cell(1,nb);      % 3D bloom coordinates
else
  ibloom = false;
  bsl = [];
  b = [];
  b3 = [];
end
%
if size(brois(1).rois(1).roi,2)<3
  itroch = false;       % No trochlear segmentations
else
  itroch = true;        % Include trochlear segmentations
end
%
f = cell(2,nrsl);       % 2D femur coordinates (1 - cartilage, 2 - bone)
t = cell(2,nrsl);       % 2D tibia coordinates (1 - cartilage, 2 - bone)
%
f3 = cell(2,nrsl);      % 3D femur coordinates (1 - cartilage, 2 - bone)
t3 = cell(2,nrsl);      % 3D tibia coordinates (1 - cartilage, 2 - bone)
%
ibone = false(nrsl,2);                 % Bone exists? (columns:  1 - femur, 2 - tibia)
maskf = false(npxt,2,nrsl);            % Mask for femoral cartilage
maskt = false(npxt,2,nrsl);            % Mask for tibial cartilage
maskb = true(npxt,nrsl);               % Mask for bloom artifacts
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
for k = 1:nrsl          % Loop through all slices
%
% Get T1 Slice Image
%
   slk = rsl(k);        % Slice number
%
   if iplt
%
     img = squeeze(v(:,:,slk));
%
     figure;
     orient landscape;
     imagesc(img);
     colormap gray;
     axis image;
     axis off;
     title({ttxt; ['Slice ' int2str(slk)]},'FontSize',16, ...
           'FontWeight','bold');
     hold on;
%
     lh = gobjects(10,1);              % Line graphic handles - maximum of five blooms per slice
     idl = false(5,1);                 % Index for legend entries
%
   end
%
% Get ROI Data for this Slice and Plot ROIs
%
   if ibloom
     nbone = 3;
   else
     nbone = 2;
   end
%
   for lb = 1:nbone     % Loop through bones - femur = 1, tibia = 2, and bloom = 3
      if lb>2
        nsurf = 1;
      else
        nsurf = 2;
      end
      for lc = 1:nsurf  % Loop through surfaces - cartilage/bloom = 1 and bone = 2
         idxl = 2*lb+lc-2;
         idxr = brois(lb).rois(lc).slice==slk;
         if any(idxr)
           if lb==1&&itroch
             nc = 3;
           elseif lb==3 % Bloom
             nc = 1;
           else
             nc = 2;
           end
           for n = 1:nc % Loop through compartments - lateral/bloom = 1, medial = 2 and trochlea = 3
              idxs = brois(lb).rois(lc).roi(n).imageno==slk;
              if any(idxs)
%
                if lb==3               % If bloom ROI
%
                  slb = find(idxs);    % Blooms on this slice
                  bn = size(slb,2);    % Number of blooms on this slice
%
                  for isl = 1:bn
%
                     idsb = slb(isl);
%
                     dat = cell2mat( ...
                                 brois(lb).rois(lc).roi(n).data(idsb)');
                     dat = fix_gap(dat);
%                    
                     dat3 = cell2mat( ...
                                brois(lb).rois(lc).roi(n).data3(idsb)');
                     dat3 = fix_gap(dat3);
%
                     b{lc,idsb} = dat;           % Bloom
                     b3{lc,idsb} = dat3;         % 3D bloom
%
                     if iplt
                       lh(isl+4) = plot(dat(:,1),dat(:,2), ...
                            lt(idxl,:),'LineWidth',0.25,'MarkerSize',2);
                       idl(5) = true;
                     end
%
% Create Bloom Artifact Mask
%
                     maskbm = ~poly2mask(dat(:,1),dat(:,2),iszs(1), ...
                                         iszs(2));         % Bloom mask
                     maskb(:,k) = maskb(:,k)&maskbm(:);    % As a vector
%
                  end   % End of for loop
%
                else    % Else not bloom ROI
%
                  dat = cell2mat(brois(lb).rois(lc).roi(n).data(idxs)');
                  dat = fix_gap(dat);
%
                  dat3 = cell2mat( ...
                                brois(lb).rois(lc).roi(n).data3(idxs)');
                  dat3 = fix_gap(dat3);
%
                  if lb==1
                    f{lc,k} = dat;     % Femur
                    f3{lc,k} = dat3;   % 3D femur
                  else
                    t{lc,k} = dat;     % Tibia
                    t3{lc,k} = dat3;   % 3D tibia
                  end
%
                  if iplt
%                    lh(idxl) = plot(dat(:,1),dat(:,2),lt(idxl,:));
                    lh(idxl) = plot(dat(:,1),dat(:,2),lt(idxl,:), ...
                                    'LineWidth',0.25,'MarkerSize',2);
                    idl(idxl) = true;
                  end
%
                end     % End of if bloom ROI
%
              end       % End of any idxs
%
           end          % End of n loop - lateral/medial/trochlear loop
         end            % End of if - any rois for this slice
      end               % End of lc loop - cartilage/bone loop
      if lb==1
        ibone(k,lb) = all(~cellfun('isempty',f(:,k)));
      end
      if lb==2
        ibone(k,lb) = all(~cellfun('isempty',t(:,k)));
      end
   end                  % End of lb loop - femur/tibia loop
%
% Add Legends and Print Slice Plots
%
   if iplt
     legend(lh(idl),legds(idl,:),'Interpreter','none');
     if iprt
       if k==1
         print('-dpsc2','-r600','-fillpage',pnam1);
       else
         print('-dpsc2','-r600','-fillpage','-append',pnam1);
       end
     end
   end
%
% Create Logical Masks for the Cartilage on this Slice
%
   if ibone(k,1)
     [maskf(:,1,k),maskf(:,2,k)] = cr_mask2f(f(:,k),iszs,dist,scal);
   end
%
   if ibone(k,2)
     [maskt(:,1,k),maskt(:,2,k)] = cr_mask2(t(:,k),iszs,dist,scal);
   end
%
% Plot ROIs
%
   if iplt
     mask1 = maskf(:,1,k)|maskt(:,1,k);     % Superficial cartilage mask
     mask1 = mask1&maskb(:,k);
     mask2 = maskf(:,2,k)|maskt(:,2,k);     % Deep cartilage mask
     mask2 = mask2&maskb(:,k);
     cmx = max(img(:));
     img1 = img-cmx-1;
     dcmx = 16*cmx/128;
     img1(mask1) = dcmx;               % Blue - Superficial
     img1(mask2) = cmx-dcmx;           % Red - Deep
%
     figure;
     orient landscape;
     imagesc(img1);
     colormap(cmap);
     caxis([-cmx cmx]);
     axis image;
     axis off;
     title({ttxt; ['Slice ' int2str(slk)]},'FontSize',16, ...
           'FontWeight','bold');
%
     if iprt
       if k==1
         print('-dpsc2','-r600','-fillpage',pnam2);
       else
         print('-dpsc2','-r600','-fillpage','-append',pnam2);
       end
     end
%
   end
%
% Close Slice Plots
%
   close all;
%
end                     % End of k loop - slice loop
%
return