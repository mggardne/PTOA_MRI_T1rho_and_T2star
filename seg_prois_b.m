%#######################################################################
%
%      * SEGmentation to PTOA Regions of Interest (ROIs) Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     CSV files to create masks for regions of interest in the femur
%     and tibial compartments.  The masks are saved in MAT files with
%     the series number and ending in "_prois.mat."
%
%     NOTES:  1.  The registered MRI MAT files must be in subject
%             directories.
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  See rd_m_dicom.m.
%
%             3.  Femoral and tibial axial segmentation CSV files for
%             the femoral and tibial axes must be in the respective
%             subdirectories.
%
%             4.  M-files InterX.m, coord_tf.m, cr_mask2.m, cr_maskfc.m,
%             cr_masktc.m, cr_mask2fc.m, cyl_fit.m, cyl_plt.m,
%             dist2cyl.m, f_cs_14.m, fem_plan.m, fix_gap.m, in_tri2d.m,
%             li_clos.m, lsect2.m, lsect2a.m, lsect3.m, lsect4.m,
%             lsect5.m, midline.m, mk_fplan_msks.m, mk_lay_msks.m,
%             mk_tplan_msks.m, mk_tri_2dc.m, mk_tri_2dfc.m,
%             mk2_tri_2d.m, mk2_tri_2df.m, mtch_ends.m, near2.m,
%             plane_fit.m, plnorm.m, plsect.m, plsect2.m, plt_datsl.m,
%             plt_masks.m, plt_maskss.m, pts2lin.m, pt2line.m,
%             rb_trnsf.m, rd_prois.m, rd_roi6hp.m, rotxyz.m, tib_plan.m,
%             tibia_cs8.m, tri2d.m, tri_area.m, and trim_seg.m must be
%             in the current directory or path.
%
%     15-Sep-2022 * Mack Gardner-Morse
%     26-Apr-2024 * Mack Gardner-Morse * Updated for baseline directory.
%

%#######################################################################
%
% Parameter for Dividing Cartilage in Half
%
% dist = 7.5;             % Maximum distance to midline in pixels
% dist = 8.0;             % Maximum distance to midline in pixels
dist = 10.0;            % Maximum distance to midline in pixels
%
% Control Plots - Femur Coordinate System? and Dividing Planes?
%
% iplt = false;            % No plot of the dividing planes
iplt = true;            % Plot the dividing planes
% iplt2 = false;          % No plot of the femur coordinate system
iplt2 = true;           % Plot the femur coordinate system
%
% Matrix to Reflect X-Axis
%
rx = eye(3);
rx(1) = -1;             % Reflect X-axis
%
% Directory with MAT Files as a String
%
dirstr = split(pwd,filesep);
%
% Get Subject
%
% subj = dirstr{end};
subj = dirstr{end-1};
%
% iskip = true;
iskip = false;
%
if ~iskip
%
% T1rho
%
id5 = 1;                % Use first spin lock time for plots
rdir = 'RHO';           % Directory for T1rho segmentations
%
% Get T1rho MAT Files in Directory
%
d = dir('T1rho_S*.mat');
mnams = {d.name}';
idr = contains(mnams,'roi','IgnoreCase',true);
mnams = mnams(~idr);
idc = contains(mnams,'chk','IgnoreCase',true);
mnams = mnams(~idc);
nmat = size(mnams,1);
%
% Loop through T1rho MAT Files
%
for m = 1:nmat
% for m = 1:1             % Right leg only
% for m = 2:nmat
%
   mnam = mnams{m};
   load(mnam,'iszs','pspcs','snt','st','v');
   v = squeeze(v(:,:,:,id5));
   fs = ['S' snt];      % Series number prefaced with a 'S'
%
% Parse Series Text for Leg
%
   if strcmpi(st(1),'L')
     leg = 'L';
     ileg = false;
     ltxt = 'Left Leg';
   else
     leg = 'R';
     ileg = true;
     ltxt = 'Right Leg';
   end
%
% Read ROIs
%
   d = dir(fullfile(rdir,[subj '_' leg '_ALL_*.csv']));
   fnam = d.name;
%
   brois = rd_prois(fullfile(rdir,fnam));   % Femur and tibia segmentations
%
   d = dir(fullfile(rdir,[subj '_' leg '_AX_FEM*.csv']));
   fnam = d.name;
%
   faxis = rd_roi6hp(fullfile(rdir,fnam));  % Femur axial segmentation
%
% Get Femur Slices
%
   rslc = brois(1).rois(1).slice;      % Femur cartilage
   rslb = brois(1).rois(2).slice;      % Femur bone
   rslbl = brois(1).rois(2).roi(1).imageno; % Femur bone lateral condyle
   rslbm = brois(1).rois(2).roi(2).imageno; % Femur bone medial condyle
   rslf = intersect(rslc,rslb);        % Ensure unique slices in sorted order
%
   nslf = size(rslf,1);                % Number of femoral slices
%
% Get Tibia Slices
%
   rslc = brois(2).rois(1).slice;      % Tibia cartilage
   rslb = brois(2).rois(2).slice;      % Tibia bone
   rslt = intersect(rslc,rslb);        % Ensure unique slices in sorted order
%
   nslt = size(rslt,1);                % Number of tibial slices
%
   idcl = brois(2).rois(1).roi(1).imageno;       % Cartilage lateral compartment
%
   idbl = brois(2).rois(2).roi(1).imageno;       % Bone lateral compartment
   idbm = brois(2).rois(2).roi(2).imageno;       % Bone medial compartment
%
   idlat = intersect(idcl,idbl)';      % Tibia lateral compartment slices
%
% Both Femur and Tibia Slices
%
   rsl = union(rslf,rslt);
   nrsl = size(rsl,1);
%
   [~,idl] = intersect(rsl,rslbl);     % Femur bone lateral condyle
   [~,idm] = intersect(rsl,rslbm);     % Femur bone medial condyle
%
   [~,idb] = intersect(rsl,rslb);      % Index for tibia bone slices
   rslb = rsl(idb);     % Get in sorted order
   [~,idbl] = intersect(rslb,idbl);    % Index for lateral tibia bone slices
   [~,idbm] = intersect(rslb,idbm);    % Index for medial tibia bone slices
%
% Create Masks that Divide Cartilage into Deep and Superficial Layers
% The first column in the logical masks is the mask for the image, the
% second column is the cartilage layers (1 - superficial, and 2 - deep),
% and the third column is the slices in RSL.
%
   [maskf,maskt,ibone,f,t,f3,t3] = mk_lay_msks(brois,rsl,nrsl,iszs, ...
                                              dist,pspcs,v,fs,subj,leg);
%
% Get Femur Coordinate System from Fitted Cylinder
%
   [xyzc,rmat] = f_cs_14(f3(2,idl)',f3(2,idm)',faxis(2).data, ...
                         faxis(1).data,iplt2);
   if iplt2
     title({'Femur Coordinate System'; ltxt},'FontSize',16, ...
           'FontWeight','bold');
   end
%
% Get Tibia Bone-Based Coordinate System and Proximal Tibia Outline
%
   d = dir(fullfile(rdir,[subj '_' leg '_AX_TIB*.csv']));
   fnam = d.name;
%
   [xyzt,rtmat,~,~,~,xyzpto] = tibia_cs8(fullfile(rdir,fnam),ileg, ...
                                         iplt2);
%
   if iplt2
     view(3);
     title({'Tibia Coordinate System'; ltxt},'FontSize',16, ...
           'FontWeight','bold');
   end
%
% Reverse X-Axis for Right Knees
%
   if ileg
     rmat = rmat*rx;    % Reflect X-axis
     rtmat = rtmat*rx;  % Reflect X-axis
   end
%
% Transform Femur Segmentations to Femur Coordinate System
%
   f3f = coord_tf(xyzc,rmat,f3(:,ibone(:,1)));
   f2 = f(:,ibone(:,1));               % Just femoral slices
%
% Transform Tibia Bone and Segmentations to Tibia Coordinate System
%
   t3tb = coord_tf(xyzt,rtmat,t3(2,idb))';
   t3t = coord_tf(xyzt,rtmat,t3(:,ibone(:,2)));
   t2 = t(:,ibone(:,2));               % Just tibial slices
%
% Create Dividing Planes for Femoral Regions of Interest (ROIs)
%
% lmp - point in lateral-medial plane
% lmv - lateral-medial plane normal vector
% ppc - point in posterior plane
% ppn - posterior plane normal vector
% ppl - point in lateral trochlea plane
% pnl - lateral trochlea plane normal vector
% ppm - point in medial trochlea plane
% pnm - medial trochlea plane normal vector
%
   f3fb = f3f(2,:)';                   % Femoral bone in femur CS
   pnam3 = [fs '_ROIs3.ps'];           % Femur dividing planes print file name
%
   [lmp,lmv,ppc,ppn,ppl,pnl,ppm,pnm] = fem_plan(f3fb,iplt,ltxt,pnam3);
%
   plan_pts = [lmp; ppc; ppl; ppm];    % Points in the dividing planes
   plan_nvs = [lmv; ppn; pnl; pnm];    % Normal vectors for the dividing planes
%
% Create Dividing Planes for Tibia Regions of Interest (ROIs)
%
% lcap - point in lateral central-anterior plane
% lcav - lateral central-anterior plane normal vector
% lpcp - point in lateral posterior-central plane
% lpcv - lateral posterior plane normal vector
% mcap - point in medial central-anterior plane
% mcav - medial central-anterior plane normal vector
% mpcp - point in medial posterior-central plane
% mpcv - medial posterior plane normal vector
%
   pnam4 = [fs '_ROIs4.ps'];           % Tibia dividing planes print file name
%
   [lcap,lcav,lpcp,lpcv,mcap,mcav,mpcp,mpcv] = tib_plan(t3tb,idbl, ...
                                                  idbm,iplt,ltxt,pnam4);
%
   tplan_pts = [lcap; lpcp; mcap; mpcp];    % Points in the dividing planes
   tplan_nvs = [lcav; lpcv; mcav; mpcv];    % Normal vectors for the dividing planes
%
   close all;
%
% Create Femoral Masks that Divide Cartilage into Regions of Interest
% (ROIs)
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
   maskfr = mk_fplan_msks(f3f,f2,rslf,nslf,plan_pts,plan_nvs,iszs, ...
                          dist,pnam3,fs,subj,leg,v);
%
%   close all;
%
% Create Tibial Masks that Divide Cartilage into Regions of Interest
% (ROIs)
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
   [masktr,t3tt] = mk_tplan_msks(t3t,t2,rslt,nslt,idlat,tplan_pts, ...
                                 tplan_nvs,iszs,dist,pnam4,fs,subj, ...
                                 leg,v);
%
   close all;
%
% Use Masks to Plot Both Femoral and Tibial ROIs
%
   pnam5 = [fs '_ROIs5.ps'];           % Tibia dividing planes print file name
%
   plt_masks(maskf,maskt,rsl,nrsl,maskfr,rslf,masktr,rslt,v, ...
             fs,subj,rdir,leg,pnam5);
%
% Save Masks, ROIS and Slice Information into MAT File
%
   savnam = [mnam(1:end-4) '_prois.mat'];
   save(savnam,'brois','f','f3','f3f','fs','ibone','leg','maskf', ...
        'maskfr','maskt','masktr','nrsl','nslf','nslt','plan_pts', ...
        'plan_nvs','rmat','rtmat','rsl','rslf','rslt','subj','t', ...
        't3','t3t','t3tt','tplan_pts','tplan_nvs','xyzc','xyzpto', ...
        'xyzt');
%
end                     % End of m loop - MAT file loop
%
end                     % End of iskip if
%
% T2star
%
rdir = 'T2S';           % Directory for T2* segmentations
%
% Get T2* MAT Files in Directory
%
d = dir('T2star_S*.mat');
mnams = {d.name}';
idr = contains(mnams,'roi','IgnoreCase',true);
mnams = mnams(~idr);
idc = contains(mnams,'chk','IgnoreCase',true);
mnams = mnams(~idc);
nmat = size(mnams,1);
%
% Loop through T2* MAT Files
%
for m = 1:nmat
% for m = 1:1
% for m = 2:nmat
%
   mnam = mnams{m};
   load(mnam,'id5','iszs','pspcs','snt','st','v');
   fs = ['S' snt];      % Series number prefaced with a 'S'
   if length(id5)>=m
     id5 = id5(m);
   else
     id5 = id5(1);
   end
   v = squeeze(v(:,:,:,id5));          % Echo time with segmentations
%
% Parse Series Text for Leg
%
   if strcmpi(st(1),'L')
     leg = 'L';
     ileg = false;
     ltxt = 'Left Leg';
   else
     leg = 'R';
     ileg = true;
     ltxt = 'Right Leg';
   end
%
% Read ROIs
%
   d = dir(fullfile(rdir,[subj '_' leg '_ALL_*.csv']));
   fnam = d.name;
%
   brois = rd_prois(fullfile(rdir,fnam));   % Femur and tibia segmentations
%
   d = dir(fullfile(rdir,[subj '_' leg '_AX_FEM*.csv']));
   fnam = d.name;
%
   faxis = rd_roi6hp(fullfile(rdir,fnam));  % Femur axial segmentation
%
% Get Femur Slices
%
   rslc = brois(1).rois(1).slice;      % Femur cartilage
   rslb = brois(1).rois(2).slice;      % Femur bone
   rslbl = brois(1).rois(2).roi(1).imageno; % Femur bone lateral condyle
   rslbm = brois(1).rois(2).roi(2).imageno; % Femur bone medial condyle
   rslf = intersect(rslc,rslb);        % Ensure unique slices in sorted order
%
   nslf = size(rslf,1);                % Number of femoral slices
%
% Get Tibia Slices
%
   rslc = brois(2).rois(1).slice;      % Tibia cartilage
   rslb = brois(2).rois(2).slice;      % Tibia bone
   rslt = intersect(rslc,rslb);        % Ensure unique slices in sorted order
%
   nslt = size(rslt,1);                % Number of tibial slices
%
   idcl = brois(2).rois(1).roi(1).imageno;       % Cartilage lateral compartment
%
   idbl = brois(2).rois(2).roi(1).imageno;       % Bone lateral compartment
   idbm = brois(2).rois(2).roi(2).imageno;       % Bone medial compartment
%
   idlat = intersect(idcl,idbl)';      % Tibia lateral compartment slices
%
% Both Femur and Tibia Slices
%
   rsl = union(rslf,rslt);
   nrsl = size(rsl,1);
%
   [~,idl] = intersect(rsl,rslbl);     % Femur bone lateral condyle
   [~,idm] = intersect(rsl,rslbm);     % Femur bone medial condyle
%
   [~,idb] = intersect(rsl,rslb);      % Index for tibia bone slices
   rslb = rsl(idb);     % Get in sorted order
   [~,idbl] = intersect(rslb,idbl);    % Index for lateral tibia bone slices
   [~,idbm] = intersect(rslb,idbm);    % Index for medial tibia bone slices
%
% Create Masks that Divide Cartilage into Deep and Superficial Layers
% The first column in the logical masks is the mask for the image, the
% second column is the cartilage layers (1 - superficial, and 2 - deep),
% and the third column is the slices in RSL.
%
   [maskf,maskt,ibone,f,t,f3,t3] = mk_lay_msks(brois,rsl,nrsl,iszs, ...
                                              dist,pspcs,v,fs,subj,leg);
%
% Get Femur Coordinate System from Fitted Cylinder
%
   [xyzc,rmat] = f_cs_14(f3(2,idl)',f3(2,idm)',faxis(2).data, ...
                         faxis(1).data,iplt2);
   if iplt2
     title({'Femur Coordinate System'; ltxt},'FontSize',16, ...
           'FontWeight','bold');
   end
%
% Get Tibia Bone-Based Coordinate System
%
   d = dir(fullfile(rdir,[subj '_' leg '_AX_TIB*.csv']));
   fnam = d.name;
%
   [xyzt,rtmat,~,~,~,xyzpto] = tibia_cs8(fullfile(rdir,fnam),ileg, ...
                                         iplt2);
%
   if iplt2
     view(3);
     title({'Tibia Coordinate System'; ltxt},'FontSize',16, ...
           'FontWeight','bold');
   end
%
% Reverse X-Axis for Right Knees
%
   if ileg
     rmat = rmat*rx;    % Reflect X-axis
     rtmat = rtmat*rx;  % Reflect X-axis
   end
%
% Transform Femur Segmentations to Femur Coordinate System
%
   f3f = coord_tf(xyzc,rmat,f3(:,ibone(:,1)));
   f2 = f(:,ibone(:,1));               % Just femoral slices
%
% Transform Tibia Bone and Segmentations to Tibia Coordinate System
%
   t3tb = coord_tf(xyzt,rtmat,t3(2,idb))';
   t3t = coord_tf(xyzt,rtmat,t3(:,ibone(:,2)));
   t2 = t(:,ibone(:,2));               % Just tibial slices
%
% Create Dividing Planes for Femoral Regions of Interest (ROIs)
%
% lmp - point in lateral-medial plane
% lmv - lateral-medial plane normal vector
% ppc - point in posterior plane
% ppn - posterior plane normal vector
% ppl - point in lateral trochlea plane
% pnl - lateral trochlea plane normal vector
% ppm - point in medial trochlea plane
% pnm - medial trochlea plane normal vector
%
   f3fb = f3f(2,:)';                   % Femoral bone in femur CS
   pnam3 = [fs '_ROIs3.ps'];           % Femur dividing planes print file name
%
   [lmp,lmv,ppc,ppn,ppl,pnl,ppm,pnm] = fem_plan(f3fb,iplt,ltxt,pnam3);
%
   plan_pts = [lmp; ppc; ppl; ppm];    % Points in the dividing planes
   plan_nvs = [lmv; ppn; pnl; pnm];    % Normal vectors for the dividing planes
%
% Create Dividing Planes for Tibia Regions of Interest (ROIs)
%
% lcap - point in lateral central-anterior plane
% lcav - lateral central-anterior plane normal vector
% lpcp - point in lateral posterior-central plane
% lpcv - lateral posterior plane normal vector
% mcap - point in medial central-anterior plane
% mcav - medial central-anterior plane normal vector
% mpcp - point in medial posterior-central plane
% mpcv - medial posterior plane normal vector
%
   pnam4 = [fs '_ROIs4.ps'];           % Tibia dividing planes print file name
%
   [lcap,lcav,lpcp,lpcv,mcap,mcav,mpcp,mpcv] = tib_plan(t3tb,idbl, ...
                                                  idbm,iplt,ltxt,pnam4);
%
   tplan_pts = [lcap; lpcp; mcap; mpcp];    % Points in the dividing planes
   tplan_nvs = [lcav; lpcv; mcav; mpcv];    % Normal vectors for the dividing planes
%
   close all;
%
% Create Femoral Masks that Divide Cartilage into Regions of Interest
% (ROIs)
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
   maskfr = mk_fplan_msks(f3f,f2,rslf,nslf,plan_pts,plan_nvs,iszs, ...
                          dist,pnam3,fs,subj,leg,v);
%
%   close all;
%
% Create Tibial Masks that Divide Cartilage into Regions of Interest
% (ROIs)
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
   [masktr,t3tt] = mk_tplan_msks(t3t,t2,rslt,nslt,idlat,tplan_pts, ...
                                 tplan_nvs,iszs,dist,pnam4,fs,subj, ...
                                 leg,v);
%
   close all;
%
% Use Masks to Plot Both Femoral and Tibial ROIs
%
   pnam5 = [fs '_ROIs5.ps'];           % Tibia dividing planes print file name
%
   plt_masks(maskf,maskt,rsl,nrsl,maskfr,rslf,masktr,rslt,v, ...
             fs,subj,rdir,leg,pnam5);
%
% Save Masks, ROIS and Slice Information into MAT File
%
   savnam = [mnam(1:end-4) '_prois.mat'];
   save(savnam,'brois','f','f3','f3f','fs','ibone','leg','maskf', ...
        'maskfr','maskt','masktr','nrsl','nslf','nslt','plan_pts', ...
        'plan_nvs','rmat','rtmat','rsl','rslf','rslt','subj','t', ...
        't3','t3t','t3tt','tplan_pts','tplan_nvs','xyzc','xyzpto', ...
        'xyzt');
%
end                     % End of m loop - MAT file loop
%
return