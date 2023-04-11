function psl_plt(v,masklay,maskroi,rsl,nrsl,rslbs,idt,tcp,nps,mxtc, ...
                 cmap,txt1,psnam)
%PSL_PLT   Plots T1rho/T2* values with the underlying images by slice 
%          within regions of interest (ROIs).
%
%          CMSL_PLT(V,MASKLAY,MASKROI,RSL,NRSL,RSLBS,IDT,TCP,NPS) Given
%          a four-dimensional matrix of T1/T2 intensities from a MRI
%          image volume, V, where the first two dimensions are an image,
%          the third dimension are the slices, and the fourth dimension
%          are the spin lock/echo times, three dimensional logical masks
%          with the first dimension being the image, the second
%          dimension being the superficial layer in the first column and
%          deep layer in the second column and the third dimension being
%          slices in a cell array of masks with the cell index to femur
%          and tibia, MASKLAY, two-dimensional logical masks with the
%          first dimension being the image, the second dimension being
%          slices in a cell array of masks with the first index to the
%          femur and tibia and the second index to lateral and medial,
%          MASKROI, array of slices with segmentations, RSL, the number
%          of slices with segmentations, NRSLS, cell array with slices
%          in the femur and tibia, RSLBS, index to the  spin lock/echo
%          time to use for plotting, IDT, cell array of T1rho/T2*
%          values, TCP, and the number of fitted pixels in each slice
%          in a cell array, NPS, plots the T1rho/T2* values with the
%          underlying images by slice within the regions of interest
%          defined by the MASKLAY and MASKROI.
%
%          PSL_PLT(V,MASKLAY,MASKROI,RSL,NRSL,RSLBS,IDT,TCP,NPS,MXTC,
%          CMAP,TXT1) Given the maximum plotting value for the color
%          scale, MXTC, a three color map, CMAP, and a text string for
%          the first line of the plot title, TXT1, plots the T1rho/T2*
%          values with a color maximum of MXTC using the color map,
%          CMAP, and using TXT1 for the first line of the plot title.
%          The default maximum value is 70.  The default color map is
%          gray for the image and jet for the T1rho/T2* values.  The
%          default first line title text is "Results Plot".
%
%          PSL_PLT(V,MASKLAY,MASKROI,RSL,NRSL,RSLBS,IDT,TCP,NPS,MXTC,
%          CMAP,TXT1,PSNAM) Given the name for a PS file, PSNAM, prints
%          the plots to the PS file.  By default, the plots are not
%          printed to a file.
%
%          NOTES:  1.  Plots the pixel output of psl_ana.m.  See
%                  psl_ana.m and mri_fitps.m.
%
%          23-Mar-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<9
  error([' *** ERROR in psl_plt:  Nine input variables are', ...
         ' required!']);
end
%
if nargin<10||isempty(mxtc)
  mxtc = 70;
end
%
if nargin<11||isempty(cmap)
%
% Default Color Map
%
  gmap = gray(128);     % Gray color map for not cartilage
  jmap = jet(128);      % Jet color map for cartilage measures
  cmap = [gmap; jmap];
end
%
if nargin<12||isempty(txt1)
  txt1 = 'Results Plot';
end
%
if nargin<13||isempty(psnam)
  isave = false;
else
  isave = true;
end
%
% Initialize Arrays
%
idxs = [2 3 2 2];        % Maximum indices for bone, compartment, ROI and layer
%
% Loop through Slices
%
for k = 1:nrsl
%
   ksl = rsl(k);        % Slice number
%
   rimg = squeeze(v(:,:,ksl,idt));     % T1/T2 data for slice and plot spin lock/echo time
%
% Scale T1rho/T2* Image to -mxtc to Zero
%      
   rimg = rimg-min(rimg(:));
   imgmx = max(rimg(:));
   rimg = mxtc*rimg./imgmx;
   rimg = rimg-(mxtc+0.01);
%
% Loop through Bones
%
   for kb = 1:2         % 1 - Femur and 2 - Tibia
%
      mskb = maskroi{kb};    % ROI masks for this bone
      rslb = rslbs{kb};      % Slices for this bone
      idxb = rslb==ksl;      % Index to slice in rslb/mskb
      if ~any(idxb)
        continue;       % No bone segmentation on this slice
      end
%
      maskb = masklay{kb};        % Layer mask for this bone
%
% Loop through Compartments
%
      for kc = 1:2      % 1 - Lateral and 2 - Medial
%
         mskc = mskb{kc};    % ROI masks for this compartment
%
% Loop through ROIs
%
         for kr = 1:3   % 1 - anterior/trochlea, 2 - central and 3 - posterior
%
            mskr = mskc{kr};
            msks = mskr(:,idxb);
%
            npps = sum(msks);
            if npps==0
              continue; % No ROI on this slice
            end
%
% Loop through Layers
%
            for kl = 1:2     % 1 - superficial, 2 - deep
%
               maskl = maskb(:,kl,k);
               mskl = msks&maskl;
%
%               npps = sum(mskl);
%
               idx = sub2ind(idxs,kl,kr,kc,kb);  % Index to T1rho/T2* results
               npsk = nps{idx};
               npsks = sum(npsk(1:k));
               npsks = (npsks-npsk(k)+1:npsks)';
%
               tcpk = tcp{idx};
%
               rimg(mskl) = tcpk(npsks);         % T1rho/T2* values
%
            end         % End of kl loop - layers loop
%
         end            % End of kr loop - ROIs loop
%
      end               % End of kc loop - compartments loop
%
   end                  % End of kb loop - bones loop
%
% Plot Slice
%
   figure;
   orient landscape;
%
   imagesc(rimg,[-mxtc mxtc]);
   colormap(cmap);
   axis image;
   axis off;
%
   if iscell(txt1)
     ttxt = {txt1{:} ['Slice ' int2str(ksl)]}';
   else
     ttxt = {txt1; ['Slice ' int2str(ksl)]};
   end
%
   title(ttxt,'FontSize',16,'FontWeight','bold');
%
   hb = colorbar;
   set(hb,'Limits',[0 mxtc]);
%
   if isave             % Print plots
     if k==1
       print('-dpsc2','-r600','-fillpage',psnam);
     else
       print('-dpsc2','-r600','-fillpage','-append',psnam);
     end
   end
%
%    close all;
%
end                     % End of k loop - slices loop
%
return