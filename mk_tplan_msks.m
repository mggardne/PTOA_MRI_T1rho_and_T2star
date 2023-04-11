function [masktr,t3tt] = mk_tplan_msks(t3t,t,rslt,nslt,idl,plan_pts, ...
                                   plan_nvs,iszs,dist,pnam,fps,sn,leg,v)
%MK_TPLAN_MSKS  Takes the PTOA knee sagittal segmentations for the tibia
%          and creates image masks for the different dividing planes 
%          between the tibia regions of interest.
%
%          MASKTR = MK_TPLAN_MSKS(T3T,T,RSLT,NSLT,IDL,PLAN_PTS,
%          PLAN_NVS,ISZS) Given a cell array with the tibia three
%          dimensional (3D) segmentations in the tibial coordinate
%          system, T3T, a cell array with the tibia two dimensional (2D)
%          segmentations in the MRI coordinate system, T, an array of
%          tibial slices with cartilage and bone segmentations, RSLT,
%          and the number of slices in RSLT, NSLT, an array of tibial
%          slices in the lateralcompartment, IDL, a matrix with the
%          coordinates of points in the dividing planes in the columns,
%          PLAN_PTS, a matrix of three dimensional vectors with the
%          vector components in the columns, PLAN_NVS, and the pixel
%          size of the MRI images, ISZS, returns the four-dimensional
%          tibial mask, MASKTR.  The first column in the logical mask
%          is the mask for the image, the second column is the dividing
%          plane (1 - lateral central-anterior, 2 - lateral
%          posterior-central, 3 - medial central-anterior, 4 - medial
%          posterior-central), the third column is the sides of plane
%          (1-anterior/2-central, 1-central/2-posterior), and the third
%          column is slices in RSLT.
%
%          [MASKTR,T2TT] = MK_TPLAN_MSKS(T3T,T,RSLT,NSLT,IDL,
%          PLAN_PTS,PLAN_NVS,ISZS) Returns a cell array with the tibia
%          three dimensional (3D) segmentations in the tibial coordinate
%          system, T3TT, with the bone segmentations truncated if they
%          extend beyond the cartilage segmentations.
%
%          MASKTR = MK_TPLAN_MSKS(T3T,T,RSLT,NSLT,IDL,PLAN_PTS,
%          PLAN_NVS,ISZS,DIST) Given a distance, DIST, ensures truncated
%          endpoints on the second boundary line are within two DIST
%          distances of the first boundary line endpoints.  See
%          mk2_tri_2d.m.
%
%          MASKTR = MK_TPLAN_MSKS(T3T,T,RSLT,NSLT,IDL,PLAN_PTS,
%          PLAN_NVS,ISZS,DIST,PNAM,FPS,SN,LEG) Given the name for a
%          Postscript (PS) file, PNAM, text for the MRI series number,
%          FPS, text for the subject number, SN, and a character for
%          the leg (L - left and R - right), LEG, plots the regions of
%          interests from the mask.  PNAM is required to generate
%          Postscript (PS) files. FPS, SN, and LEG are optional texts
%          that are used to generate plot titles.  PNAM must be input
%          (may be empty to not save to a file) to generate plots.
%
%          MASKTR = MK_TPLAN_MSKS(T3T,T,RSLT,NSLT,IDL,PLAN_PTS,
%          PLAN_NVS,ISZS,DIST,PNAM,FPS,SN,LEG,V) Given the image matrix
%          (third dimension is the slices), V, generates two-dimensional
%          (2D) plots of the femur ROIs overlaid on the image slices.
%
%          NOTES:  1.  PNAM must be input (may be empty to not save to
%                  a file) to generate plots.
%
%                  2.  PNAM is required to generate a Postscript (PS)
%                  file.
%
%                  3.  FPS, SN, and LEG are used to title the plots.
%                  Any underscores in FPS are replaced with spaces.
%
%                  4.  The M-files cr_maskt.m, cr_masktc.m, decomp.m,
%                  in_tri2d.m, lsect2.m, lsect2a.m, mk2_tri_2d.m,
%                  mk2_tri_2dc.m, mtch_ends.m, near2.m, plane_fit.m,
%                  rb_trnsf.m, and trim_seg.m the current directory or
%                  path.
%
%          18-Jan-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check Number of Inputs
%
if (nargin<8)
  error([' *** ERROR in MK_TPLAN_MSKS:  At least eight (8) inputs', ...
         ' are required!']);
end
%
if nargin<9||isempty(dist)
  dist = Inf;           % No distance checking
end
%
if (nargin>9)
  iplt = true;
else
  iplt = false;
end
%
if (nargin>=10)&&~isempty(pnam)
  iprt = true;
else
  iprt = false;
end
%
if iprt
%
  ttxt = strrep(fps,'_',' ');
%  
  if exist('sn','var')&&~isempty(sn)
    ttxt = [ttxt ', Subject ' sn];
  end
%  
  if exist('leg','var')&&~isempty(leg)
    if strcmpi(leg,'L')
      ttxt = [ttxt ', Left Leg'];
    else
      ttxt = [ttxt ', Right Leg'];
    end
  end
end
%
if nargin==14&&~isempty(v)
  iplt2 = true;
else
  iplt2 = false;
end
%
% Check Inputs
%
if ~iscell(t3t)
  error([' *** ERROR in MK_TPLAN_MSKS:  First input must be a', ...
         ' cell array with the tibial segmentations!']);
end
%
% Initialize ROIs Mask and Truncated Slice Coordinates
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
npx = prod(iszs);       % Number of pixels
nplanes = size(plan_pts,1);            % Number of planes
%
masktr = false(npx,nplanes,2,nslt);
%
t3tt = t3t;             % Truncated 3D slice coordinates
%
% Setup Color Maps, and Legend Strings
%
if iplt2
%
  gmap = gray(128);     % Gray color map for not cartilage
  jmap = jet(128);      % Jet color map for cartilage
  cmap = [gmap; jmap];
  cscal = [0.16 0.33 0.5 0.6 0.7 0.9];
%
end
%
% Setup Patch Colors and Figures for 3D Plots
%
if iplt
%
  colrs = 'mbcgry';     % Patch colors
%  legds = {'Lateral Anterior'; 'Lateral Central'; ...
%           'Lateral Posterior'; 'Medial Anterior'; ...
%           'Medial Central'; 'Medial Posterior'};
%
  hf = gobjects(nplanes,1);
%
% 3D Plane Plot Figures
%
  plnams = {'Lateral Central-Anterior'; 'Lateral Posterior-Central'; ...
            'Medial Central-Anterior';  'Medial Posterior-Central'};
%
  for k = 1:nplanes
     hf(k) = figure;
     orient landscape;
     hold on;
  end
%
end
%
%  Loop through the Tibial Slices
%
for k = 1:nslt
%
   ksl = rslt(k);       % Slice number
   ilat = any(idl==ksl);               % True for lateral slices
%
% Truncate Bone Segmentations based on the Lengths of the Cartilage
% Segmentations
%
   t3tt(:,k) = trim_seg(t3t(:,k),dist);
%
% Get Cartilage and Bone Segmentations
%
   t3tc = t3tt{1,k};    % Cartilage on this slice
   t3tb = t3tt{2,k};    % Bone on this slice
   tc = t{1,k};         % 2D cartilage on this slice
   tb = t{2,k};         % 2D bone on this slice
%
% Get Transformation from 3D Femoral Coordinates to 2D Pixel Coordinates
%
% s - scale
% r - rotation matrix
% t - translation vector
%
   nptsb = size(tb,1);
   [s,r,tv] = decomp(t3t{2,k},[tb zeros(nptsb,1)]);
%
% Loop through the Tibial Planes
%
   if ilat
     iplane1 = 1;       % First lateral plane
     iplane2 = 2;       % Last lateral plane
   else
     iplane1 = 3;       % First medial plane
     iplane2 = 4;       % Last medial plane
   end
%
   for l = iplane1:iplane2
%
% Get Cartilage and Bone Intersections with the Plane
%
      plan_pt = plan_pts(l,:);         % Point for this plane
      plan_nv = plan_nvs(l,:);         % Normal vector for this plane
%
      [ptc,~,idxc] = plsect2(plan_pt,plan_nv,t3tc);   % Cartilage
      [ptb,~,idxb] = plsect2(plan_pt,plan_nv,t3tb);   % Bone
%
% If No Intersection - Determine Side of Plane
% If an Intersection - Split Slice between Sides of Plane
%
      if isempty(ptc)   % No cartilage intersection
%
        if ~isempty(ptb)               % Bone only intersection
%
          warning([' *** Warning in MK_TPLAN_MSKS:  Tibia', ...
                   ' cartilage and bone both do not have', ...
                   ' intersections on slice ' int2str(ksl) '!']);
%
% Truncate Bone
%
          nb = size(t3tb,1);
          if idxb>nb/2
            t3tbt = [t3tb(1:idxb(1),:); ptb];
          else
            t3tbt = [ptb; t3tb(idxb(1)+1:nb,:)];
          end
%
% Transform Back to Pixel Coordinates
%
          nt3d = size(t3tbt,1);        % Number of points
          tbt = s*t3tbt*r+repmat(tv,nt3d,1);     % Transform to pixel coordinates
          tbt = tbt(:,1:2);
%
% Create Mask
%
          side = mean(t3tc)*plan_nv';
          planex = plan_pt*plan_nv';
          if side<planex               % Central/Posterior
            masktr(:,l,2,k) = cr_maskt({tc; tbt},iszs,dist);
            p = 2;
          else                         % Anterior/Central
            masktr(:,l,1,k) = cr_maskt({tc; tbt},iszs,dist);
            p = 1;
          end
%
% Plot 3D Regions of Interest (ROIs) for this Plane and this Slice
%
          if iplt
%
            figure(hf(l));
%
            dr1(2,:) = t3tbt(end,:)-t3tbt(1,:);
            dr1(1,:) = t3tc(end,:)-t3tc(1,:);
            pdot = dr1(1,:)*dr1(2,:)';
            if pdot>0
              t3tbt = flipud(t3tbt);
            end
            pxyz = [t3tc; t3tbt];
            patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),colrs(2*p), ...
                  'EdgeColor','none','FaceAlpha',0.33);
%
          end
%
        else            % No intersections (neither cartilage nor bone)
%
          side = mean(t3tc)*plan_nv';
          planex = plan_pt*plan_nv';
          if side<planex               % Central/Posterior
            masktr(:,l,2,k) = cr_maskt({tc; tb},iszs,dist);
            p = 2;
          else                         % Anterior/Central
            masktr(:,l,1,k) = cr_maskt({tc; tb},iszs,dist);
            p = 1;
          end
%
% Plot 3D Regions of Interest (ROIs) for this Plane and this Slice
%
          if iplt
%
            figure(hf(l));
%
            dr1(2,:) = t3tb(end,:)-t3tb(1,:);
            dr1(1,:) = t3tc(end,:)-t3tc(1,:);
            pdot = dr1(1,:)*dr1(2,:)';
            if pdot>0
              t3tb = flipud(t3tb);
            end
            pxyz = [t3tc; t3tb];
            patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),colrs(2*p), ...
                  'EdgeColor','none','FaceAlpha',0.33);
%
          end

        end
%
      else              % Cartilage intersection
%
        if isempty(ptb)                % Cartilage only intersection
%
          warning([' *** Warning in MK_TPLAN_MSKS:  Tibia', ...
                   ' cartilage and bone both do not have', ...
                   ' intersections on slice ' int2str(ksl) '!']);
%
% Truncate Cartilage
%
          nc = size(t3tc,1);
          if idxc>nc/2
            t3tct = [t3tc(1:idxc(1),:); ptc];
          else
            t3tct = [ptc; t3tc(idxc(1)+1:nc,:)];
          end
%
% Transform Back to Pixel Coordinates
%
          nt3d = size(t3tct,1);        % Number of points
          tct = s*t3tct*r+repmat(tv,nt3d,1);     % Transform to pixel coordinates
          tct = tct(:,1:2);
%
% Create Mask
%
          side = mean(t3tct)*plan_nv';
          planex = plan_pt*plan_nv';
          if side<planex               % Central/Posterior
            masktr(:,l,2,k) = cr_maskt({tct; tb},iszs,dist);
            p = 2;
          else                         % Anterior/Central
            masktr(:,l,1,k) = cr_maskt({tct; tb},iszs,dist);
            p = 1;
          end
%
% Plot 3D Regions of Interest (ROIs) for this Plane and this Slice
%
          if iplt
%
            figure(hf(l));
%
            dr1(2,:) = t3tb(end,:)-t3tb(1,:);
            dr1(1,:) = t3tct(end,:)-t3tct(1,:);
            pdot = dr1(1,:)*dr1(2,:)';
            if pdot>0
              t3tb = flipud(t3tb);
            end
            pxyz = [t3tct; t3tb];
            patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),colrs(2*p), ...
                  'EdgeColor','none','FaceAlpha',0.33);
%
          end
%
        else            % Cartilage and bone intersections 
%
          t3tc1 = [t3tc(1:idxc(1),:); ptc];
          t3tc2 = [ptc; t3tc(idxc(1)+1:end,:)];
%
          side1 = mean(t3tc1)*plan_nv';
          side2 = mean(t3tc2)*plan_nv';
          if side1<side2
            idsc(1,2) = 1;
            idsc(1,1) = 2;
          else
            idsc(1,2) = 2;
            idsc(1,1) = 1;
          end
%
          t3tb1 = [t3tb(1:idxb(1),:); ptb];
          t3tb2 = [ptb; t3tb(idxb(1)+1:end,:)];
%
          side1 = mean(t3tb1)*plan_nv';
          side2 = mean(t3tb2)*plan_nv';
          if side1<side2
            idsb(1,2) = 1;
            idsb(1,1) = 2;
          else
            idsb(1,2) = 2;
            idsb(1,1) = 1;
          end
%
% Transform Back to Pixel Coordinates
%
          [tr,tr3,dr] = rb_trnsf(s,r,tv,t3tc1,t3tb1,t3tc2,t3tb2);
%
% Find Any Overlapping Points in the Bone Segmentation
%
          [~,mtchb] = mtch_ends(tr);
          clip = ~[(mtchb(1:2)|mtchb(3:4))'; ...
                 (mtchb([1 3])|mtchb([2 4]))'];  % Clip not matching ends
%
% Create Masks
%
          for p = 1:2   % Two plane directions
%
             idxc = idsc==p;
             idxb = idsb==p;
             pdat = {tr{1,idxc}; tr{2,idxb}};
             clp = clip(idxb,:);
             masktr(:,l,p,k) = cr_masktc(pdat,iszs,clp,dist);   % Create mask
%
% Plot 3D Regions of Interest (ROIs) for this Plane and this Slice
%
             if iplt
               figure(hf(l));
               pdat3 = {tr3{1,idxc}; tr3{2,idxb}};
               pdot = dr(1,:,idxc)*dr(2,:,idxb)';
               if pdot>0
                 pdat3{1} = flipud(pdat3{1});
               end
               pxyz = cell2mat(pdat3);
               patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),colrs(2*p), ...
                     'EdgeColor','none','FaceAlpha',0.33);
             end
%
          end           % End of p loop - plane directions
%
        end             % End of cartilage intersection
%
      end               % End of if intersection
%
   end                  % End of l loop - planes
%
% Plot 2D Regions of Interest (ROIs) for this Slice
%
   if iplt2
%
% Get ROI Masks
%
     maskla = masktr(:,1,1,k);         % Lateral anterior
     masklc = masktr(:,1,2,k)&masktr(:,2,1,k);   % Lateral central
     masklp = masktr(:,2,2,k);         % Lateral posterior
%
     maskma = masktr(:,3,1,k);         % Medial anterior
     maskmc = masktr(:,3,2,k)&masktr(:,4,1,k);   % Medial central
     maskmp = masktr(:,4,2,k);         % Medial posterior
%
% Get MRI Image with ROIs
%
     img = squeeze(v(:,:,ksl));
     cmx = max(img(:));
     img = img-cmx-1;
     cscalx = cscal*cmx;
     img(maskla) = cscalx(1);          % Blue - Lateral anterior
     img(masklc) = cscalx(2);          % Cyan - Lateral central
     img(masklp) = cscalx(3);          % Green - Lateral posterior
     img(maskma) = cscalx(6);          % Red - Medial anterior
     img(maskmc) = cscalx(5);          % Orange - Medial central
     img(maskmp) = cscalx(4);          % Yellow - Medial posterior
%
% Plot 2D ROIs
%
     hfp = figure;
     orient landscape;
     imagesc(img);
     colormap(cmap);
     caxis([-cmx cmx]);
     axis image;
     axis off;
     title({ttxt; ['Slice ' int2str(ksl)]},'FontSize',16, ...
           'FontWeight','bold');
%
     if iprt
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
%
     close(hfp);        % Close 2D pixel slice figure
%
   end                  % End of if iplt2
%
end                     % End of k loop - slices
%
% Finish and Print 3D Figures
%
if iplt
  for k = 1:nplanes
     figure(hf(k));
     view(3);
     axis equal;
     title({ttxt; [plnams{k} ' Plane']},'FontSize',16, ...
           'FontWeight','bold');
     axlim = axis;
     axlim = reshape(axlim,2,3);
     plan_pt = plan_pts(k,:);          % Point for this plane
     plan_nv = plan_nvs(k,:);          % Normal vector for this plane
     [xp,yp,zp] = plane_plt(plan_pt,plan_nv,axlim,2);
     surf(xp,yp,zp,'EdgeColor','r','FaceColor',[0.7 0.7 0.7], ...
          'FaceAlpha',0.33);
     view(3);
     if iprt
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
%
  end
end
%
return