function plt_masks(maskf,maskt,rsl,nrsl,maskfr,rslf,masktr,rslt,v, ...
                   fps,sn,rdir,leg,pnam)
%PLT_MASKS Plots the cartilage regions of interests for the PTOA knee
%          sagittal segmentations for the femur and tibia.
%
%          PLT_MASKS(MASKF,MASKT,RSL,NRSL,MASKFR,RSLF,MASKTR,RSLT,V)
%          Given three-dimensional logical mask arrays for the deep and
%          superficial cartilage layers for the femur, MASKF, and tibia,
%          MASKT, a list of image slice numbers for the masks, RSL, the
%          number of slices, NRSL, four-dimensional logical mask array
%          for the regions of interests (ROIs) in the femur, MASKFR, a
%          list of image slice numbers for the femur ROIs mask, RSLF,
%          four-dimensional logical mask array for the ROIs in the
%          tibia, MASKTR, a list of image slice numbers for the tibia
%          ROIs mask, RSLT, the image matrix (third dimension is the
%          slices), V, generates plots of the slices with the ROIs and
%          layers color coded.
%
%          PLT_MASKS(MASKF,MASKT,RSL,NRSL,MASKFR,RSLF,MASKTR,RSLT,V,
%          FPS,SN,RDIR,LEG) Given the text for the MRI series number,
%          FPS, text for the subject number, SN, the analysis
%          directory, RDIR, and a character for the leg (L - left and
%          R - right), LEG, generates plot titles.
%
%          PLT_MASKS(MASKF,MASKT,RSL,NRSL,MASKFR,RSLF,MASKTR,RSLT,V,
%          FPS,SN,RDIR,LEG,PNAM) Given the name for a Postscript (PS)
%          file, PNAM, prints the plots to the file PNAM.
%
%          NOTES:  1.  PNAM is required to generate a Postscript (PS)
%                  file.
%
%                  2.  FPS, SN, and LEG are optional and are used to
%                  title the plots. Any underscores in FPS are replaced
%                  with spaces.
%
%          27-Jan-2022 * Mack Gardner-Morse 
%

%#######################################################################
%
% Check Number of Inputs
%
if (nargin<9)
  error([' *** ERROR in PLT_MASKS:  At least nine (9) inputs', ...
         ' are required!']);
end
%
if (nargin>=10)&&~isempty(fps)
%
  ttxt = strrep(fps,'_',' ');
%  
  if exist('sn','var')&&~isempty(sn)
    ttxt = [ttxt ', Subject ' sn];
  end
%  
  if exist('rdir','var')&&~isempty(rdir)
    if strcmpi(rdir,'RHO')
      ttxt = [ttxt ', T1\rho'];
    else
      ttxt = [ttxt ', T2*'];
    end
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
if (nargin==14)&&~isempty(pnam)
  iprt = true;
else
  iprt = false;
end
%
% Check Inputs
%
if ~islogical(maskf)
  error([' *** ERROR in PLT_MASKS:  First input must be a', ...
         ' logical array!']);
end
%
% Setup Color Maps for 2D Plots
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
cscal = [0.1562 0.3125 0.5078 0.6055 0.7031 0.8984;
         0.0781 0.3516 0.5469 0.6641 0.7422 0.9766];
%
%  Loop through the MRI Slices with Regions of Interest (ROIs)
%
for k = 1:nrsl
%
   ksl = rsl(k);        % Slice number
   idxf = rslf==ksl;    % Index to slice in rslf/maskfr
   idxt = rslt==ksl;    % Index to slice in rslt/masktr
%
% Get MRI Image with ROIs
%
   img = squeeze(v(:,:,ksl));
   cmx = max(img(:));
   img = img-cmx-1;
   cscalx = cscal*cmx;
%
% Get Femur ROI Masks
%
   if any(idxf)
%
     masklt = maskfr(:,1,1,idxf)&maskfr(:,3,1,idxf);
     masklc = maskfr(:,1,1,idxf)&maskfr(:,2,1,idxf)&maskfr(:,3,2,idxf);
     masklp = maskfr(:,1,1,idxf)&maskfr(:,2,2,idxf);
%
     maskmt = maskfr(:,1,2,idxf)&maskfr(:,4,1,idxf);
     maskmc = maskfr(:,1,2,idxf)&maskfr(:,2,1,idxf)&maskfr(:,4,2,idxf);
     maskmp = maskfr(:,1,2,idxf)&maskfr(:,2,2,idxf);
%
% Loop through the Cartilage Layers (1 - Superficial and 2 - Deep)
%
     for l = 1:2
%
        maskltl = masklt&maskf(:,l,k);
        masklcl = masklc&maskf(:,l,k);
        masklpl = masklp&maskf(:,l,k);
        maskmtl = maskmt&maskf(:,l,k);
        maskmcl = maskmc&maskf(:,l,k);
        maskmpl = maskmp&maskf(:,l,k);
%
        m = 3-l;
        img(maskltl) = cscalx(m,1);    % Blue/Dark Blue - Lateral trochlea
        img(masklcl) = cscalx(m,2);    % Light Blue/Cyan - Lateral central
        img(masklpl) = cscalx(m,3);    % Green/Lime - Lateral posterior
        img(maskmtl) = cscalx(m,6);    % Red/Brown - Medial trochlea
        img(maskmcl) = cscalx(m,5);    % Dark Tan/Orange - Medial central
        img(maskmpl) = cscalx(m,4);    % Yellow/Light Tan - Medial posterior
%
     end
   end
%
% Get Tibia ROI Masks
%
   if any(idxt)
%
     maskla = masktr(:,1,1,idxt);      % Lateral anterior
     masklc = masktr(:,1,2,idxt)&masktr(:,2,1,idxt);  % Lateral central
     masklp = masktr(:,2,2,idxt);      % Lateral posterior
%
     maskma = masktr(:,3,1,idxt);      % Medial anterior
     maskmc = masktr(:,3,2,idxt)&masktr(:,4,1,idxt);  % Medial central
     maskmp = masktr(:,4,2,idxt);      % Medial posterior
%
% Loop through the Cartilage Layers (1 - Superficial and 2 - Deep)
%
     for l = 1:2
%
        masklal = maskla&maskt(:,l,k);
        masklcl = masklc&maskt(:,l,k);
        masklpl = masklp&maskt(:,l,k);
        maskmal = maskma&maskt(:,l,k);
        maskmcl = maskmc&maskt(:,l,k);
        maskmpl = maskmp&maskt(:,l,k);
%
        img(masklal) = cscalx(l,1);    % Dark Blue/Blue - Lateral anterior
        img(masklcl) = cscalx(l,2);    % Cyan/Light Blue - Lateral central
        img(masklpl) = cscalx(l,3);    % Lime/Green - Lateral posterior
        img(maskmal) = cscalx(l,6);    % Brown/Red - Medial anterior
        img(maskmcl) = cscalx(l,5);    % Orange/Dark Tan - Medial central
        img(maskmpl) = cscalx(l,4);    % Light Tan/Yellow - Medial posterior
%
     end
   end
%
% Plot ROIs
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
     if k==1
       print('-dpsc2','-r600','-fillpage',pnam);
     else
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
   end
%
   close(hfp);          % Close 2D pixel slice figure
%
end                     % End of k loop - slices
%
return