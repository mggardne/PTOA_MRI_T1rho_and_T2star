function pcmprt_plt(v,masklay,maskroi,rsls,nrsls,rsl,rslbs,idt,tcp, ...
                    nps,mxtc,cmap,txt1,psnam)
%PCMPRT_PLT Plots T1rho/T2* values with the underlying images by slice 
%          within regions of interest (ROIs).
%
%          CMPRT_PLT(V,MASKLAY,MASKROI,RSLS,NRSLS,RSL,RSLBS,IDT,TCP,NPS)
%          Given a four-dimensional matrix of T1/T2 intensities from a MRI image
%          volume, V, where the first two dimensions are an image, the
%          third dimension are the slices, and the fourth dimension are
%          the spin lock/echo times, three dimensional logical masks
%          with the first dimension being the image, the second
%          dimension being the superficial layer in the first column and
%          deep layer in the second column and the third dimension being
%          slices in a cell array of masks with the first index to the
%          lateral and medial compartments and the second index to the
%          femur and tibia, MASK, a cell array with the slices within
%          each compartment, RSLS, the number of slices in each
%          compartment, NRSLS, index to the  spin lock/echo time to use
%          for plotting, IDT, cell array of T1rho/T2* values, TCP, and
%          the number of fitted pixels in each slice within the
%          compartments in a cell array, NPS, plots the T1rho/T2*
%          values with the underlying images by slice within the
%          regions of interest defined by the MASK.
%
%          CMPRT_PLT(V,MASK,RSLS,NRSLS,IDT,TCP,NPS,MXTC,CMAP,TXT1) Given
%          the maximum plotting value for the color scale, MXTC, a three
%          color map, CMAP, and a text string for the first line of the
%          plot title, TXT1, plots the T1rho/T2* values with a color
%          maximum of MXTC using the color map, CMAP, and using TXT1 for
%          the first line of the plot title.  The default maximum value
%          is 70.  The default color map is gray for the image and jet
%          for the T1rho/T2* values.  The default first line title text
%          is "Results Plot".
%
%          CMPRT_PLT(V,MASK,RSLS,NRSLS,IDT,TCP,NPS,MXTC,CMAP,TXT1,PSNAM)
%          Given the name for a PS file, PSNAM, prints the plots to the
%          PS file.  By default, the plots are not printed.
%
%          NOTES:  1.  Plots the pixel output of cmprt_ana.m.  See
%                  cmprt_ana.m and mri_fitr2.m.
%
%          16-Jun-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<7
  error([' *** ERROR in pcmprt_plt:  Seven input variables are', ...
         ' required!']);
end
%
if nargin<8||isempty(mxtc)
  mxtc = 70;
end
%
if nargin<9||isempty(cmap)
%
% Default Color Map
%
  gmap = gray(128);     % Gray color map for not cartilage
  jmap = jet(128);      % Jet color map for cartilage measures
  cmap = [gmap; jmap];
end
%
if nargin<10||isempty(txt1)
  txt1 = 'Results Plot';
end
%
if nargin<11||isempty(psnam)
  isave = false;
else
  isave = true;
end
%
% Compartment Labels
%
tcmprts = {'Lateral'; 'Medial'};
%
% Initialize Arrays
%
idxs = [3 2 2 2];       % Maximum indices for compartment, bone, layer, and ROI
%
% Loop through Compartments
%
for kc = 1:2
%
   mskc = maskroi{kc};  % ROI masks for this compartment
   rslr = rsls{kc};     % Slices for this compartment
   nrslr = nrsls(kc);   % Number of slices in this compartment
   tcmprt = [tcmprts{kc} ' Compartment'];
%
% Loop through Slices
%
   for ks = 1:nrslr
%
      slk = rslr(ks);   % Slice
%
% Get Slice Image
%
      rimg = squeeze(v(:,:,slk,idt));  % T1/T2 data for slice and plot spin lock/echo time
%
% Scale T1rho/T2* Image to -mxtc to Zero
%      
      rimg = rimg-min(rimg(:));
      imgmx = max(rimg);
      rimg = mxtc*rimg./imgmx;
      rimg = rimg-(mxtc+0.01);
%
% Loop through Bone
%
      for kb = 1:2
%
         maskb = masklay{kb};          % Layer mask for this bone
         mskb = mskc{kb};              % ROI mask for this bone
%
         rslb = rslbs{kb};             % Slices for this bone
%
         idxl = rsl==slk;              % Index to layer slices in rsl/masklay
         idxb = rslb==slk;             % Index to bone slices in rslbs/maskroi
%
% Loop through Layer
%
         for kl = 1:2   % 1 - superficial and 2 - deep
%
            maskl = maskb(:,kl,idxl);
%
            for kr = 1:3               % 1 - anterior/trochlea, 2 - central and 3 - posterior
%
               mskr = mskb{kr};
               msks = mskr(:,idxb);
               mskl = msks&maskl;
%
               idx = sub2ind(idxs,kr,kl,kb,kc);  % Index to T1rho/T2* results
               npsk = nps{idx};
               npsks = sum(npsk(1:ks));
               npsks = (npsks-npsk(ks)+1:npsks)';
%
               tcpk = tcp{idx};
%
               rimg(mskl) = tcpk(npsks);         % T1rho/T2* values
%
            end         % End of kr loop - ROIs loop
%
         end            % End of kl loop - layers loop
      end               % End of kb loop - bones loop
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
      title({txt1; ['Slice ' int2str(slk)]; tcmprt},'FontSize',16, ...
            'FontWeight','bold');
%
      hb = colorbar;
      set(hb,'Limits',[0 mxtc]);
%
      if isave          % Print plots
        if kr==1&&ks==1
          print('-dpsc2','-r600','-fillpage',psnam);
        else
          print('-dpsc2','-r600','-fillpage','-append',psnam);
        end
      end
%
   end                  % End of ks loop - slices loop
%
%    close all;
%
end                     % End of kr loop - compartments loop
%
return