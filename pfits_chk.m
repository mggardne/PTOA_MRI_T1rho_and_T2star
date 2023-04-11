%#######################################################################
%
%                     * PTOA FITS CHecK Program *
%
%          M-File which reads subject 003 right leg T1rho (Series 1401)
%     and left leg T2* (Series 501) and plots the T1rho or T2* data.
%     For the femur, only the deep layer of the trochlea and posterior
%     regions, and the superficial layer of the center region are
%     plotted.  For the tibia, only the superficial layer of the
%     anterior and posterior regions, and the deep layer of the center
%     region.  The posterior regions are outlined in red pixels and the
%     lateral center regions are outlined in blue pixels.
%
%          Plots are output to Postscript files.
%
%     NOTES:  1.  MAT files from dicom_chk3.m for Series 501 and 1401
%             must be in the current directory or path.
%
%             2.  MAT files from mri_fits.m must be in the directory
%             \..\Results\mri_fitps\.
%
%     03-Apr-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Map
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
%
% Parameters
%
trmx = 100;             % Maximum valid T1rho result
trmn = 0;               % Minimum valid T1rho result
tsmx = 100;             % Maximum valid T2* result
tsmn = 0;               % Minimum valid T2* result
%
mxtr = 80;              % Maximum scale on T1rho plots
mxts = 75;              % Maximum scale on T2* plots
%
% Get Structured Element for Eroding the Image by One Pixel
%
sec = strel(logical([0 1 0; 1 1 1; 0 1 0]));     % Crossed lines
%
% MS-Excel Output Spreadsheet
%
ifirst = true;          % First write to file
xlsnams = 'pfits_chk_S';     % Start of results spreadsheet file name
xlstyp = '.xlsx';
hdrs1 = {'Subject' 'Result' 'Leg' 'Bone' 'Comprt' 'ROI' 'Layer', ...
         'Slice'};
hdrs2 = {'Pixels' 'ValidPixChk' 'ValidPix' 'MeanCheck' 'Mean' ...
         'MeanDiff' 'MinDiff' 'MaxDiff' 'RSM Diff'};
%
% PS Output File
%
psnams = 'pfits_chk_S';      % Start of PS file name
pstyp = '.ps';               % PS file type
%
% Load T1rho Data
%
% iskip = true;
iskip = false;
if ~iskip
load('T1rho_S1401.mat','snt','v');
load('T1rho_S1401_prois.mat','maskf','maskfr','maskt','masktr', ... 
     'rsl','rslf','rslt');
load T1rho_S1401_chk;
s = load(fullfile('..','Results','mri_fitps','mri_fitps.mat'), ...
         't1r_nps','t1r_respx');
nps = s.t1r_nps;
tcp = s.t1r_respx;
%
% Include Series Number in Output File Names
%
xlsnam = [xlsnams snt xlstyp];
psnam = [psnams snt pstyp];
%
idt = 1;                % Spin lock time for plots - 1 = 0 ms spin lock time
%
% Get Femur ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
fmasklt = squeeze(maskfr(:,1,1,:)&maskfr(:,3,1,:));   % Lateral trochlea
fmasklc = squeeze(maskfr(:,1,1,:)&maskfr(:,2,1,:)& ...
                 maskfr(:,3,2,:));     % Lateral central
fmasklp = squeeze(maskfr(:,1,1,:)&maskfr(:,2,2,:));   % Lateral posterior
%
fmaskl = {fmasklt; fmasklc; fmasklp};  % Lateral femur
%
fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));   % Medial trochlea
fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                 maskfr(:,4,2,:));     % Medial central
fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));   % Medial posterior
%
fmaskm = {fmaskmt; fmaskmc; fmaskmp};  % Medial femur
%
% Get Tibia ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
tmaskla = squeeze(masktr(:,1,1,:));    % Lateral anterior
tmasklc = squeeze(masktr(:,1,2,:)&masktr(:,2,1,:));   % Lateral central
tmasklp = squeeze(masktr(:,2,2,:));    % Lateral posterior
%
tmaskl = {tmaskla; tmasklc; tmasklp};  % Lateral tibia
%
tmaskma = squeeze(masktr(:,3,1,:));    % Medial anterior
tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));   % Medial central
tmaskmp = squeeze(masktr(:,4,2,:));    % Medial posterior
%
tmaskm = {tmaskma; tmaskmc; tmaskmp};  % Medial tibia
%
% Combine Masks into Cell Arrays
%
masklay = {maskf; maskt};    % Combine femur and tibia masks
%
fmask = {fmaskl; fmaskm};    % Combine femur masks
tmask = {tmaskl; tmaskm};    % Combine tibia masks
%
maskroi = {fmask; tmask};    % Combine femur and tibia masks  
%
rslbs = {rslf; rslt};        % Combine femur and tibia slices
%
% Loop through Slices
%
ksubj = 1;              % Subject index
kres = 0;               % Result identifier (0 - T1rho, and 1 - T2*)
kleg = 2;               % Leg index (1 - left, and 2 - right)
%
nsl = size(plt_sl,1);
%
for p = 1:nsl
%
   ks = plt_sl(p);
   k = find(rsl==ks);   % Index to layer masks
%
% Scale T1rho/T2* Image to -mxtr to Zero
%      
   img = squeeze(v(:,:,ks,idt));  % T1 data for slice
%
   img = img-min(img(:));
   imgmx = max(img(:));
   img = mxtr*img./imgmx;
   img = img-(mxtr+0.01);
%
%    rimg1 = img;         % Full slice check values
   rimg2 = img;         % Masked ROI values
%
% Get Check Values
%
   T1rhonlsk = squeeze(T1rhonls(:,:,p));
   T1rhonlsk = T1rhonlsk(:);
%
% Loop through Bones
%
   for kb = 1:2         % 1 - Femur and 2 - Tibia
%
      mskb = maskroi{kb};    % ROI masks for this bone
      rslb = rslbs{kb};      % Slices for this bone
      idxb = rslb==ks;       % Index to slice in rslb/mskb
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
% mod(kr+kb,2)+1 produces alternating deep and superficial layers in
% the femur and tibia:
%
%         femur  tibia  ROIs
%  kr/kb    1      2
%            layers 1 - deep, 2 - superficial
%    1      1      2    trochlea/anterior
%    2      2      1    center
%    3      1      2    posterior
%
            for kl = mod(kr+kb,2)+1    % 1 - superficial, 2 - deep
%
               maskl = maskb(:,3-kl,k);
               mskl = msks&maskl;
%
               T1rhonlsl = T1rhonlsk(mskl);
               tcpk = tcp{ksubj,kleg,kb,kc,kr,kl};
%
               npsk = nps{ksubj,kleg,kb,kc,kr,kl};
               npsks = sum(npsk(1:p));
               npsks = (npsks-npsk(p)+1:npsks)';
%
               tcpkl = tcpk(npsks);
%
%                rimg1(mskl) = T1rhonlsl;     % T1rho/T2* values
               rimg2(mskl) = tcpkl;         % T1rho/T2* values
%
               if kr==3
                 mske = mskl;
                 msk = imerode(mskl,sec);
                 mske(msk) = false;    % Mask of edge
                 rimg2(mske) = mxts;   % Red edge for posterior
               end
%
               if kr==2&&kc==1
                 mske = mskl;
                 msk = imerode(mskl,sec);
                 mske(msk) = false;    % Mask of edge
                 rimg2(mske) = 0;      % Blue edge for lateral center
               end
%
% Comparison of T1rho/T2* Values
%
               npx = npsk(p);
%
               idv1 = T1rhonlsl>=trmn&T1rhonlsl<=trmx;
               idv2 = tcpkl>=trmn&tcpkl<=trmx;
               iddv = idv1&idv2;
%
               npxv1 = sum(idv1);
               npxv2 = sum(idv2);
%
               vdiff = tcpkl-T1rhonlsl;
               vdiff = vdiff(iddv);    % Compare only valid differences
%
               t1ra = mean(T1rhonlsl(idv1));
               tcpa = mean(tcpkl(idv2));
%
               vda = mean(vdiff);
               vdmn = min(vdiff);
               vdmx = max(vdiff);
               vdrms = sqrt((vdiff'*vdiff)/sum(iddv));
%
% Combine Identifiers
%
               id = [kb kc kr kl]-1;
               ids = [ksubj kres kleg-1 id ks];
%
% Create and Write Table of Results
%
               t1 = array2table(ids,'VariableNames',hdrs1);
               t2 = table(npx,npxv1,npxv2,t1ra,tcpa,vda,vdmn,vdmx, ...
                          vdrms,'VariableNames',hdrs2);
               t = [t1 t2];
%
               if ifirst
                 writetable(t,xlsnam,'WriteMode','replacefile');
                 ifirst = false;
               else
                 writetable(t,xlsnam,'WriteMode','append', ...
                            'WriteVariableNames',false);
               end
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
   orient tall;
%
%    subplot(2,1,1);
%    imagesc(rimg1,[-mxtr mxtr]);
%    colormap(cmap);
%    axis image;
%    axis off;
%
%    ttxt = {['Subject 003, Series ' snt]; ['T1\rho, Right Leg,', ...
%            ' Slice ' int2str(ks)]; 'Full Slice Check Values'};
   ttxt = {['Subject 003, Series ' snt]; ['T1\rho, Right Leg,', ...
           ' Slice ' int2str(ks)]};
%
%    title(ttxt,'FontSize',12,'FontWeight','bold');
%
%    hb = colorbar;
%    set(hb,'Limits',[0 mxtr]);
%
%    subplot(2,1,2);
   imagesc(rimg2,[-mxtr mxtr]);
   colormap(cmap);
   axis image;
   axis off;
%
%     ttxt = 'Masked ROI Values';
%
   title(ttxt,'FontSize',12,'FontWeight','bold');
%
   hb = colorbar;
   set(hb,'Limits',[0 mxtr]);
%
   if p==1
     print('-dpsc2','-r600','-fillpage',psnam);
   else
     print('-dpsc2','-r600','-fillpage','-append',psnam);
   end
%
   close all;
%
end                     % End of p loop - slices loop
end                     % End of if iskip
%
% Load T2* Data
%
load('T2star_S501.mat','snt','v');
load('T2star_S501_prois.mat','maskf','maskfr','maskt','masktr', ... 
     'rsl','rslf','rslt');
load T2star_S501_1001_chk;
s = load(fullfile('..','Results','mri_fitps','mri_fitps.mat'), ...
         't2s_nps','t2s_respx');
nps = s.t2s_nps;
tcp = s.t2s_respx;
%
% Include Series Number in Output File Names
%
ifirst = true;          % First write to file
xlsnam = [xlsnams snt xlstyp];
psnam = [psnams snt pstyp];
%
idt = 3;                % Echo time for plots - 3 = 5 ms echo time
%
% Get Femur ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
fmasklt = squeeze(maskfr(:,1,1,:)&maskfr(:,3,1,:));   % Lateral trochlea
fmasklc = squeeze(maskfr(:,1,1,:)&maskfr(:,2,1,:)& ...
                 maskfr(:,3,2,:));     % Lateral central
fmasklp = squeeze(maskfr(:,1,1,:)&maskfr(:,2,2,:));   % Lateral posterior
%
fmaskl = {fmasklt; fmasklc; fmasklp};  % Lateral femur
%
fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));   % Medial trochlea
fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                 maskfr(:,4,2,:));     % Medial central
fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));   % Medial posterior
%
fmaskm = {fmaskmt; fmaskmc; fmaskmp};  % Medial femur
%
% Get Tibia ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
tmaskla = squeeze(masktr(:,1,1,:));    % Lateral anterior
tmasklc = squeeze(masktr(:,1,2,:)&masktr(:,2,1,:));   % Lateral central
tmasklp = squeeze(masktr(:,2,2,:));    % Lateral posterior
%
tmaskl = {tmaskla; tmasklc; tmasklp};  % Lateral tibia
%
tmaskma = squeeze(masktr(:,3,1,:));    % Medial anterior
tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));   % Medial central
tmaskmp = squeeze(masktr(:,4,2,:));    % Medial posterior
%
tmaskm = {tmaskma; tmaskmc; tmaskmp};  % Medial tibia
%
% Combine Masks into Cell Arrays
%
masklay = {maskf; maskt};    % Combine femur and tibia masks
%
fmask = {fmaskl; fmaskm};    % Combine femur masks
tmask = {tmaskl; tmaskm};    % Combine tibia masks
%
maskroi = {fmask; tmask};    % Combine femur and tibia masks  
%
rslbs = {rslf; rslt};        % Combine femur and tibia slices
%
% Loop through Slices
%
ksubj = 1;              % Subject index
kres = 1;               % Result identifier (0 - T1rho, and 1 - T2*)
kleg = 1;               % Leg index (1 - left, and 2 - right)
%
nsl = size(plt_sl,1);
%
for p = 1:nsl
%
   ks = plt_sl(p);
   k = find(rsl==ks);   % Index to layer masks
%
% Scale T1rho/T2* Image to -mxts to Zero
%      
   img = squeeze(v(:,:,ks,idt));  % T1 data for slice
%
   img = img-min(img(:));
   imgmx = max(img(:));
   img = mxts*img./imgmx;
   img = img-(mxts+0.01);
%
%    rimg1 = img;         % Full slice check values
   rimg2 = img;         % Masked ROI values
%
% Get Check Values
%
   T2starnlsk = squeeze(T2starnls(:,:,p));
   T2starnlsk = T2starnlsk(:);
%
% Loop through Bones
%
   for kb = 1:2         % 1 - Femur and 2 - Tibia
%
      mskb = maskroi{kb};    % ROI masks for this bone
      rslb = rslbs{kb};      % Slices for this bone
      idxb = rslb==ks;       % Index to slice in rslb/mskb
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
% mod(kr+kb,2)+1 produces alternating deep and superficial layers in
% the femur and tibia:
%
%         femur  tibia  ROIs
%  kr/kb    1      2
%            layers 1 - deep, 2 - superficial
%    1      1      2    trochlea/anterior
%    2      2      1    center
%    3      1      2    posterior
%
            for kl = mod(kr+kb,2)+1    % 1 - superficial, 2 - deep
%
               maskl = maskb(:,3-kl,k);
               mskl = msks&maskl;
%
               if sum(mskl)==0
                 continue; % No ROI on this slice
               end
%
               T2starnlsl = T2starnlsk(mskl);
               tcpk = tcp{ksubj,kleg,kb,kc,kr,kl};
%
               npsk = nps{ksubj,kleg,kb,kc,kr,kl};
               npsks = sum(npsk(1:p));
               npsks = (npsks-npsk(p)+1:npsks)';
%
               tcpkl = tcpk(npsks);
%
%                rimg1(mskl) = T2starnlsl;    % T1rho/T2* values
               rimg2(mskl) = tcpkl;         % T1rho/T2* values
%
               if kr==3
                 mske = mskl;
                 msk = imerode(mskl,sec);
                 mske(msk) = false;    % Mask of edge
                 rimg2(mske) = mxts;   % Red edge for posterior
               end
%
               if kr==2&&kc==1
                 mske = mskl;
                 msk = imerode(mskl,sec);
                 mske(msk) = false;    % Mask of edge
                 rimg2(mske) = 0;      % Blue edge for lateral center
               end
%
% Comparison of T1rho/T2* Values
%
               npx = npsk(p);
%
               idv1 = T2starnlsl>=trmn&T2starnlsl<=trmx;
               idv2 = tcpkl>=trmn&tcpkl<=trmx;
               iddv = idv1&idv2;
%
               npxv1 = sum(idv1);
               npxv2 = sum(idv2);
%
               vdiff = tcpkl-T2starnlsl;
               vdiff = vdiff(iddv);    % Compare only valid differences
%
               t1ra = mean(T2starnlsl(idv1));
               tcpa = mean(tcpkl(idv2));
%
               vda = mean(vdiff);
               vdmn = min(vdiff);
               vdmx = max(vdiff);
               vdrms = sqrt((vdiff'*vdiff)/sum(iddv));
%
% Combine Identifiers
%
               id = [kb kc kr kl]-1;
               ids = [ksubj kres kleg-1 id ks];
%
% Create and Write Table of Results
%
               t1 = array2table(ids,'VariableNames',hdrs1);
               t2 = table(npx,npxv1,npxv2,t1ra,tcpa,vda,vdmn,vdmx, ...
                          vdrms,'VariableNames',hdrs2);
               t = [t1 t2];
%
               if ifirst
                 writetable(t,xlsnam,'WriteMode','replacefile');
                 ifirst = false;
               else
                 writetable(t,xlsnam,'WriteMode','append', ...
                            'WriteVariableNames',false);
               end
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
%    orient tall;
   orient landscape;
%
%    subplot(2,1,1);
%    imagesc(rimg1,[-mxts mxts]);
%    colormap(cmap);
%    axis image;
%    axis off;
%
%    ttxt = {['Subject 003, Series ' snt]; ['T2*, Left Leg,', ...
%            ' Slice ' int2str(ks)]; 'Full Slice Check Values'};
   ttxt = {['Subject 003, Series ' snt]; ['T2*, Left Leg,', ...
           ' Slice ' int2str(ks)]};
%
%    title(ttxt,'FontSize',12,'FontWeight','bold');
%
%    hb = colorbar;
%    set(hb,'Limits',[0 mxts]);
%
%    subplot(2,1,2);
   imagesc(rimg2,[-mxts mxts]);
   colormap(cmap);
   axis image;
   axis off;
%
%    ttxt = 'Masked ROI Values';
%
   title(ttxt,'FontSize',12,'FontWeight','bold');
%
   hb = colorbar;
   set(hb,'Limits',[0 mxts]);
%
   if p==1
     print('-dpsc2','-r600','-fillpage',psnam);
   else
     print('-dpsc2','-r600','-fillpage','-append',psnam);
   end
%
   close all;
%
end                     % End of p loop - slices loop
