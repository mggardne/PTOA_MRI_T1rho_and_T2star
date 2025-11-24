function maskfr = mk_fplan_msks(f3f,f,rslf,nslf,plan_pts,plan_nvs, ...
                                iszs,dist,pnam,fps,sn,leg,v)
%MK_FPLAN_MSKS  Takes the PTOA knee sagittal segmentations for the femur
%          and creates image masks for the different dividing planes 
%          between the femur regions of interest.
%
%          MASKFR = MK_FPLAN_MSKS(F3F,F,RSLF,NSLF,PLAN_PTS,PLAN_NVS,
%          ISZS) Given a cell array with the femur three dimensional
%          (3D) segmentations in the femoral coordinate system, F3F, a
%          cell array with the femur two dimensional (2D) segmentations
%          in the MRI coordinate system, F, an array of femoral slices
%          with cartilage and bone segmentations, RSLF, and the number
%          of slices in RSLF, NSLF, a matrix with the coordinates of
%          points in the dividing planes in the columns, PLAN_PTS, a
%          matrix of three dimensional vectors with the vector
%          components in the columns, PLAN_NVS, and the pixel size of
%          the MRI images, ISZS, returns the four-dimensional femoral
%          mask, MASKFR.  The first column in the logical mask is the
%          mask for the image, the second column is the dividing plane
%          (1 - lateral-medial, 2 - posterior-center, 3 - lateral
%          trochlea, and 4 - medial trochlea), the third column is the
%          sides of plane (1-lateral/2-medial, 1-center/2-posterior,
%          1-lateral trochlea/2-center, and 1-medial trochlea/2-center),
%          and the third column is slices in RSLF.
%
%          MASKFR = MK_FPLAN_MSKS(F3F,F,RSLF,NSLF,PLAN_PTS,PLAN_NVS,
%          ISZS,DIST) Given a distance, DIST, ensures truncated
%          endpoints on the second boundary line are within two DIST
%          distances of the first boundary line endpoints.  See
%          mk2_tri_2df.m.
%
%          MASKFR = MK_FPLAN_MSKS(F3F,F,RSLF,NSLF,PLAN_PTS,PLAN_NVS,
%          ISZS,DIST,PNAM,FPS,SN,LEG) Given the name for a Postscript
%          (PS) file, PNAM, text for the MRI series number, FPS, text
%          for the subject number, SN, and a character for the leg
%          (L - left and R - right), LEG, plots the regions of
%          interests from the mask.  PNAM is required to generate a
%          Postscript (PS) file. FPS, SN, and LEG are optional texts
%          that are used to generate plot titles.  PNAM must be input
%          (may be empty to not save to a file) to generate plots.
%
%          MASKFR = MK_FPLAN_MSKS(F3F,F,RSLF,NSLF,PLAN_PTS,PLAN_NVS,
%          ISZS,DIST,PNAM,FPS,SN,LEG,V) Given the image matrix (third
%          dimension is the slices), V, generates two-dimensional (2D)
%          plots of the femur ROIs overlaid on the image slices.
%
%          NOTES:  1.  PNAM must be input (may be empty to not save to
%                  a file) to generate plots.
%
%                  2.  PNAM is required to a generate Postscript (PS)
%                  file.
%
%                  3.  FPS, SN, and LEG are used to title the plots.
%                  Any underscores in FPS are replaced with spaces.
%
%                  4.  The M-files cr_maskf.m, cr_maskfc.m, decomp.m,
%                  in_tri2d.m, lsect2.m, lsect2a.m, mk2_tri_2df.m,
%                  mk2_tri_2dfc.m, mtch_ends.m, near2.m, plsect2.m, and
%                  rb_trnsf.m must be in the current directory or path.
%
%          09-Jan-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check Number of Inputs
%
if (nargin<7)
  error([' *** ERROR in MK_FPLAN_MSKS:  At least seven (7) inputs', ...
         ' are required!']);
end
%
if nargin<8||isempty(dist)
  dist = Inf;           % No distance checking
end
%
if (nargin>8)
  iplt = true;
else
  iplt = false;
end
%
if (nargin>=9)&&~isempty(pnam)
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
if nargin==13&&~isempty(v)
  iplt2 = true;
else
  iplt2 = false;
end
%
% Check Inputs
%
if ~iscell(f3f)
  error([' *** ERROR in MK_FPLAN_MSKS:  First input must be a', ...
         ' cell array with the femoral segmentations!']);
end
%
% Initialize ROIs Mask
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
npx = prod(iszs);       % Number of pixels
nplanes = size(plan_pts,1);            % Number of planes
%
maskfr = false(npx,nplanes,2,nslf);
%
% Setup Color Maps for 2D Plots
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
%  legds = {'Lateral Trochlea'; 'Lateral Central'; ...
%           'Lateral Posterior'; 'Medial Trochlea'; ...
%           'Medial Central'; 'Medial Posterior'};
%
  hf = gobjects(nplanes,1);
%
% 3D Plane Plot Figures
%
  plnams = {'Lateral-Medial'; 'Posterior-Center'; ... 
            'Lateral Trochlea'; 'Medial Trochlea'};
%
  for k = 1:nplanes
     hf(k) = figure;
     orient landscape;
     hold on;
  end
%
end
%
%  Loop through the Femoral Slices
%
for k = 1:nslf
%
   ksl = rslf(k);       % Slice number
   idx = rslf==ksl;     % Index to slice in rslf/f/f3f/mask
%
% Get Cartilage and Bone Segmentations
%
   f3fc = f3f{1,idx};               % Cartilage on this slice
   f3fb = f3f{2,idx};               % Bone on this slice
   fc = f{1,idx};                   % 2D cartilage on this slice
   fb = f{2,idx};                   % 2D bone on this slice
%
% Get Transformation from 3D Femoral Coordinates to 2D Pixel Coordinates
%
% s - scale
% r - rotation matrix
% t - translation vector
%
   nptsb = size(fb,1);
   [s,r,t] = decomp(f3fb,[fb zeros(nptsb,1)]);
%
% Loop through the Femoral Planes
%
   for l = 1:nplanes
%
% Get Cartilage and Bone Intersections with the Plane
%
      plan_pt = plan_pts(l,:);         % Point for this plane
      plan_nv = plan_nvs(l,:);         % Normal vector for this plane
%
      [ptc,~,idxc] = plsect2(plan_pt,plan_nv,f3fc);   % Cartilage
      [ptb,~,idxb] = plsect2(plan_pt,plan_nv,f3fb);   % Bone
%
% If No Intersection - Determine Side of Plane
% If an Intersection - Split Slice between Sides of Plane
%
      if isempty(ptc)   % No cartilage intersection
%
        if ~isempty(ptb)               % Bone only intersection
%
          warning([' *** Warning in MK_FPLAN_MSKS:  Femur', ...
                   ' cartilage and bone both do not have', ...
                   ' intersections on slice ' int2str(ksl) '!']);
%
% Truncate Bone
%
          nb = size(f3fb,1);
          if idxb>nb/2
            f3fbt = [f3fb(1:idxb(1),:); ptb];
          else
            f3fbt = [ptb; f3fb(idxb(1)+1:nb,:)];
          end
%
% Transform Back to Pixel Coordinates
%
          nf3d = size(f3fbt,1);        % Number of points
          fbt = s*f3fbt*r+repmat(t,nf3d,1);      % Transform to pixel coordinates
          fbt = fbt(:,1:2);
%
% Create Mask
%
          side = sign(mean(f3fbt)*plan_nv');
          if side<0     % Medial/Posterior/Center
            maskfr(:,l,2,k) = cr_maskf({fc; fbt},iszs,dist);
            p = 2;
          else          % Lateral/Center/Trochlea
            maskfr(:,l,1,k) = cr_maskf({fc; fbt},iszs,dist);
            p = 1;
          end
%
% Plot 3D Regions of Interest (ROIs) for this Plane and this Slice
%
          if iplt
%
            figure(hf(l));
%
            dr1(2,:) = f3fbt(end,:)-f3fbt(1,:);
            dr1(1,:) = f3fc(end,:)-f3fc(1,:);
            pdot = dr1(1,:)*dr1(2,:)';
            if pdot>0
              f3fbt = flipud(f3fbt);
            end
            pxyz = [f3fc; f3fbt];
            patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),colrs(2*p), ...
                  'EdgeColor','none','FaceAlpha',0.33);
%
          end
%
        else            % No intersections (neither cartilage nor bone)
%
          side = sign(mean(f3fc)*plan_nv');
          if side<0     % Medial/Posterior/Center
            maskfr(:,l,2,k) = cr_maskf({fc; fb},iszs,dist);
            p = 2;
          else          % Lateral/Center/Trochlea
            maskfr(:,l,1,k) = cr_maskf({fc; fb},iszs,dist);
            p = 1;
          end
%
% Plot 3D Regions of Interest (ROIs) for this Plane and this Slice
%
          if iplt
%
            figure(hf(l));
%
            dr1(2,:) = f3fb(end,:)-f3fb(1,:);
            dr1(1,:) = f3fc(end,:)-f3fc(1,:);
            pdot = dr1(1,:)*dr1(2,:)';
            if pdot>0
              f3fb = flipud(f3fb);
            end
            pxyz = [f3fc; f3fb];
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
          warning([' *** Warning in MK_FPLAN_MSKS:  Femur', ...
                   ' cartilage and bone both do not have', ...
                   ' intersections on slice ' int2str(ksl) '!']);
%
% Truncate Cartilage
%
          nc = size(f3fc,1);
          if idxc>nc/2
            f3fct = [f3fc(1:idxc(1),:); ptc];
          else
            f3fct = [ptc; f3fc(idxc(1)+1:nc,:)];
          end
%
% Transform Back to Pixel Coordinates
%
          nf3d = size(f3fct,1);        % Number of points
          fct = s*f3fct*r+repmat(t,nf3d,1);      % Transform to pixel coordinates
          fct = fct(:,1:2);
%
% Create Mask
%
          side = sign(mean(f3fct)*plan_nv');
          if side<0     % Medial/Posterior/Center
            maskfr(:,l,2,k) = cr_maskf({fct; fb},iszs,dist);
            p = 2;
          else          % Lateral/Center/Trochlea
            maskfr(:,l,1,k) = cr_maskf({fct; fb},iszs,dist);
            p = 1;
          end
%
% Plot 3D Regions of Interest (ROIs) for this Plane and this Slice
%
          if iplt
%
            figure(hf(l));
%
            dr1(2,:) = f3fb(end,:)-f3fb(1,:);
            dr1(1,:) = f3fct(end,:)-f3fct(1,:);
            pdot = dr1(1,:)*dr1(2,:)';
            if pdot>0
              f3fb = flipud(f3fb);
            end
            pxyz = [f3fct; f3fb];
            patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),colrs(2*p), ...
                  'EdgeColor','none','FaceAlpha',0.33);
%
          end
%
        else            % Cartilage and bone intersections 
%
          f3fc1 = [f3fc(1:idxc(1),:); ptc];
          f3fc2 = [ptc; f3fc(idxc(1)+1:end,:)];
%
          side = sign(mean(f3fc1)*plan_nv');
          if side<0     % Medial/Posterior/Center
            idsc(1,2) = 1;
            idsc(1,1) = 2;
          else          % Lateral/Center/Trochlea
            idsc(1,2) = 2;
            idsc(1,1) = 1;
          end
%
          f3fb1 = [f3fb(1:idxb(1),:); ptb];
          f3fb2 = [ptb; f3fb(idxb(1)+1:end,:)];
%
          side = sign(mean(f3fb1)*plan_nv');
          if side<0     % Medial/Posterior/Center
            idsb(1,2) = 1;
            idsb(1,1) = 2;
          else          % Lateral/Center/Trochlea
            idsb(1,2) = 2;
            idsb(1,1) = 1;
          end
%
% Transform Back to Pixel Coordinates
%
          [fr,fr3,dr] = rb_trnsf(s,r,t,f3fc1,f3fb1,f3fc2,f3fb2);
%
% Find Any Overlapping Points in the Bone Segmentation
%
          [~,mtchb] = mtch_ends(fr);
          clip = ~[(mtchb(1:2)|mtchb(3:4))'; ...
                 (mtchb([1 3])|mtchb([2 4]))'];  % Clip not matching ends
%
% Create Masks
%
          for p = 1:2   % Two plane directions
%
             idxc = idsc==p;
             idxb = idsb==p;
             pdat = {fr{1,idxc}; fr{2,idxb}};
             clp = clip(idxb,:);
             try
                maskfr(:,l,p,k) = cr_maskfc(pdat,iszs,clp,dist);     % Create mask
             catch
                warning([' *** Warning in MK_FPLAN_MSKS:  Femur', ...
                   ' cartilage and bone do not form an ROI', ...
                   ' on slice ' int2str(ksl) '!']);
%                keyboard
             end
%
% Plot 3D Regions of Interest (ROIs) for this Plane and this Slice
%
             if iplt
               figure(hf(l));
               pdat3 = {fr3{1,idxc}; fr3{2,idxb}};
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
     masklt = maskfr(:,1,1,k)&maskfr(:,3,1,k);
     masklc = maskfr(:,1,1,k)&maskfr(:,2,1,k)&maskfr(:,3,2,k);
     masklp = maskfr(:,1,1,k)&maskfr(:,2,2,k);
%
     maskmt = maskfr(:,1,2,k)&maskfr(:,4,1,k);
     maskmc = maskfr(:,1,2,k)&maskfr(:,2,1,k)&maskfr(:,4,2,k);
     maskmp = maskfr(:,1,2,k)&maskfr(:,2,2,k);
%
% Get MRI Image with ROIs
%
     img = squeeze(v(:,:,ksl));
     cmx = max(img(:));
     img = img-cmx-1;
     cscalx = cscal*cmx;
     img(masklt) = cscalx(1);          % Blue - Lateral trochlea
     img(masklc) = cscalx(2);          % Cyan - Lateral central
     img(masklp) = cscalx(3);          % Green - Lateral posterior
     img(maskmt) = cscalx(6);          % Red - Medial trochlea
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
     if iprt
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
%
  end
end
%
return