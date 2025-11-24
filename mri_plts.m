%#######################################################################
%
%                      * MRI PLoT Slice Program *
%
%          M-File which reads the segmentation MAT files and results
%     MAT file to plot the T2* results overlaid on the MRI data for
%     particular slices.  The results are saved as PNG image files.
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "0".
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "prois".  See rd_m_dicom.m and
%             seg_prois.m.
%
%             3.  M-file psl_plts.m must be in the current directory or
%             path.
%
%             4.  Only for T2* plots.  Currently setup to get T2* plots
%             for female lateral compartment slices for T1rho/T2* paper.
%
%     27-Sep-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Map
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
%
% Initialize Parameters
%
mxtr = 80;              % Maximum scale on T1rho plots
mxts = 75;              % Maximum scale on T2* plots
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','mri_fitps');        % Results directory
%
pngnam = fullfile(resdir,'mri_fitps_');     % Start of PNG file name
% pngtyp = '.png';        % PNG file type
%
% Get T2* Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Bone - 1 = femur and 2 = tibia
%   Index 4 - Compartment - 1 = lateral and 2 = medial
%   Index 5 - ROI - 1 = anterior/trochlea, 2 = central and 3 = posterior
%   Index 6 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and compartment analysis variables are:
%        1 = superficial and 2 = deep.
%
load(fullfile(resdir,'mri_fitps.mat'),'t2s_res','t2s_npx', ...
     't2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
%
% Get Subject Directories
%
sdirs = dir('0*');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
%
% Get Index to First 20 Subjects without Subject 019
%
id19 = ~contains(sdirs,'019');
sdirs = sdirs(id19);
%
% Subjects for Plots
%
idsubj = contains(sdirs,pattern({'018';'031';'051'}));     % Subjects 18, 31, and 51
sdirs = sdirs(idsubj);
nsubj = size(sdirs,1);
%
idsubj = find(idsubj);
%
% Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Bone - 1 = femur and 2 = tibia
%   Index 4 - Compartment - 1 = lateral and 2 = medial
%   Index 5 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 6 - Layer - 1 = deep and 2 = superficial
%
sls = [16 33];          % Left and right leg lateral compartment slices
%
% Loop through Subjects
%
for ks = 1:nsubj
% for ks = nsubj:nsubj
% for ks = 1:5
   kss = idsubj(ks);    % Index to subject results
%
% Get Subject Directory (Name) and Number
%
   sdir = sdirs{ks};    % Current subject directory (name)
   subj = eval(sdir);   % Subject number
   subjtxt = ['Subject ' sdir];
%
   pngnams = [pngnam sdir];            % Add subject to PNG file name
%
% Get T1rho MAT Files in Directory
%
   ido = false;         % Skip T1rho
%    ido = true;          % Do T1rho
%
   if ido
%
     d = dir(fullfile(sdir,'T1rho_S*.mat'));
     roinams = {d.name}';
     idr = contains(roinams,'proi','IgnoreCase',true);     % Masks
%
     rhonams = roinams(~idr);
     idx = contains(rhonams,'chk','IgnoreCase',true);      % Check files
     rhonams = rhonams(~idx);          % Image MAT files
     nrho = size(rhonams,1);
%
     roinams = roinams(idr);           % ROI MAT files
     nroi = size(roinams,1);
%
     if nrho~=nroi
       error([' *** ERROR in mri_fitps:  Number of T1rho MAT files', ...
              ' does not match the number of ROI MAT files!']);
     end
     clear nroi;
%
% T1rho Identifier
%
     ires = 0;          % ires = 0 - T1rho, ires = 1 - T2*
     restxt = [subjtxt ' T1\rho'];
     idt = 1;           % Spin lock/echo time for plots - 1 = 0 ms spin lock time
%
     psnamr = [pngnams '_T1R_'];       % Add result type to PS file name
%
% Loop through T1rho MAT Files
%
     for km = 1:nrho
%
% Load Data
%
        rhonam = rhonams{km};
        load(fullfile(sdir,rhonam),'iszs','nslt','scmx','sns', ...
             'snt','splt','st','v');
        npix = prod(iszs);   % Number of pixels in an image
        fs = ['S' snt];      % Series number prefaced with a 'S'
%
        idm = contains(roinams,rhonam(1:end-4)); % Get matching file
        roinam = roinams{idm};
        load(fullfile(sdir,roinam),'maskf','maskfr', ...
             'maskt','masktr','nrsl','rsl','rslf','rslt');
%
% Parse Series Text for Leg
%
        if strcmpi(st(1),'L')
          leg = 'L';
          legtxt = ' Left Leg';
          ileg = 0;     % Coding for leg
        else
          leg = 'R';
          legtxt = ' Right Leg';
          ileg = 1;     % Coding for leg
        end
%
% Add Leg to PS File Name
%
        psnamf = [psnamr leg pstyp];   % Add leg to PS file name
%
% Get Femur ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
        fmasklt = squeeze(maskfr(:,1,1,:)&maskfr(:,3,1,:));     % Lateral trochlea
        fmasklc = squeeze(maskfr(:,1,1,:)&maskfr(:,2,1,:)& ...
                         maskfr(:,3,2,:));       % Lateral central
        fmasklp = squeeze(maskfr(:,1,1,:)&maskfr(:,2,2,:));     % Lateral posterior
%
        fmaskl = {fmasklt; fmasklc; fmasklp};    % Lateral femur
%
        fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));     % Medial trochlea
        fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                         maskfr(:,4,2,:));       % Medial central
        fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));     % Medial posterior
%
        fmaskm = {fmaskmt; fmaskmc; fmaskmp};    % Medial femur
%
% Get Tibia ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
        tmaskla = squeeze(masktr(:,1,1,:));      % Lateral anterior
        tmasklc = squeeze(masktr(:,1,2,:)&masktr(:,2,1,:));     % Lateral central
        tmasklp = squeeze(masktr(:,2,2,:));      % Lateral posterior
%
        tmaskl = {tmaskla; tmasklc; tmasklp};    % Lateral tibia
%
        tmaskma = squeeze(masktr(:,3,1,:));      % Medial anterior
        tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));     % Medial central
        tmaskmp = squeeze(masktr(:,4,2,:));      % Medial posterior
%
        tmaskm = {tmaskma; tmaskmc; tmaskmp};    % Medial tibia
%
% Combine Masks into Cell Arrays
%
        masklay = {maskf; maskt};      % Combine femur and tibia masks
%
        fmask = {fmaskl; fmaskm};      % Combine femur masks
        tmask = {tmaskl; tmaskm};      % Combine tibia masks
        maskroi = {fmask; tmask};      % Combine femur and tibia masks
%
        rslbs = {rslf; rslt};          % Combine femur and tibia slices
%
% Do Slice Analysis
%
%         [tc,~,rss,npx,id,tcp,~,rssp,nps] = psl_ana(v,masklay, ...
%                      maskroi,rsl,nrsl,rslbs,splt,nslt,fun,init,tr0,opt);
%         na = size(tc,1);               % Number of results
%
% Get Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Bone - 1 = femur and 2 = tibia
%   Index 4 - Compartment - 1 = lateral and 2 = medial
%   Index 5 - ROI - 1 = anterior/trochlea, 2 = central and 3 = posterior
%   Index 6 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and analysis variables are:
%        1 = superficial and 2 = deep.
%
        for ka = 1:na
           t1r_res(ks,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = tc(ka);
           t1r_npx(ks,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = npx(ka);
           t1r_rss(ks,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                   id(ka,4)+1) = rss(ka);
           t1r_respx{ks,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = tcp{ka};
           t1r_rsspx{ks,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = rssp{ka};
           t1r_nps{ks,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                     id(ka,4)+1} = nps{ka};
        end
%
% Plot Results
%
        sid = [restxt legtxt];    % Subject, result (T1rho/T2*), and leg
%         psl_plt(v,masklay,maskroi,rsl,nrsl,rslbs,idt,tcp,nps,mxtr, ...
%                 cmap,sid,psnamf);
%
% Get Statistics on Pixel Results
%
        npxv = zeros(na,1);            % Number of valid results
        tcpm = zeros(na,1);            % Mean
        tcpmn = zeros(na,1);           % Minimum
        tcpmx = zeros(na,1);           % Maximum
        tcpsd = zeros(na,1);           % SD
%
        for ka = 1:na
           idv = tcp{ka}>=trmn&tcp{ka}<=trmx;
           npxv(ka) = sum(idv);        % Number of valid results
           if npxv(ka)>0
             tcpv = tcp{ka}(idv);      % Valid T1rho values
             tcpm(ka) = mean(tcpv);    % Mean
             tcpmn(ka) = min(tcpv);    % Minimum
             tcpmx(ka) = max(tcpv);    % Maximum
             tcpsd(ka) = std(tcpv);    % SD
           end
        end
%
        tcpcov = 100*tcpsd./tcpm;      % Coefficient of variation
        tcpcov(isnan(tcpcov)) = 0;     % Catch any NaNs
%
% Combine Identifiers
%
        ids = [subj ires ileg];        % MAT file identifiers
        ids = repmat(ids,na,1);
        ids = [ids id];                % All identifiers
%
% Create and Write Table of Results
%
        t1 = array2table(ids,'VariableNames',hdrs1);
        t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                   'VariableNames',hdrs2);
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
     end                % End of km loop - T1rho MAT file loop
%
     close all;         % Close all plot windows
%
   end                  % End of ido - Skip T1rho?
%
% Get T2* MAT Files in Directory
%
   d = dir(fullfile(sdir,'T2star_S*.mat'));
   roinams = {d.name}';
   idr = contains(roinams,'proi','IgnoreCase',true);       % Masks
%
   starnams = roinams(~idr);
   idx = contains(starnams,'chk','IgnoreCase',true);       % Check files
   starnams = starnams(~idx);          % Image MAT files
   nstar = size(starnams,1);
%
   roinams = roinams(idr);             % ROI MAT files
   nroi = size(roinams,1);
%
   if nstar~=nroi
     error([' *** ERROR in mri_fitps:  Number of T2* MAT files', ...
            ' does not match the number of ROI MAT files!']);
   end
   clear nroi;
%
% T2* Identifier
%
   restxt = [subjtxt ' T2*'];
   idt = 3;             % Spin lock/echo time for plots - 3 = 5 ms echo time
%
   pngnamr = [pngnams '_T2S_'];        % Add result type to PNG file name
%
% Loop through T2* MAT Files
%
   for km = 1:nstar
%
% Load Data
%
      starnam = starnams{km};
      load(fullfile(sdir,starnam),'scmx','sns','snt','st','v');
      fs = ['S' snt];        % Series number prefaced with a 'S'
%
      idm = contains(roinams,starnam(1:end-4));    % Get matching file
      roinam = roinams{idm};
      load(fullfile(sdir,roinam),'maskf','maskfr', ...
           'maskt','masktr','nrsl','rsl','rslf','rslt');
%
% Parse Series Text for Leg
%
      if strcmpi(st(1),'L')
        leg = 'L';
        legtxt = ' Left Leg';
        ileg = 0;       % Coding for leg
      else
        leg = 'R';
        legtxt = ' Right Leg';
        ileg = 1;       % Coding for leg
      end
%
% Add Leg to PNG File Name
%
      pngnamf = [pngnamr leg];         % Add leg to PNG file name
%
% Get Femur ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
      fmasklt = squeeze(maskfr(:,1,1,:)&maskfr(:,3,1,:));  % Lateral trochlea
      fmasklc = squeeze(maskfr(:,1,1,:)&maskfr(:,2,1,:)& ...
                       maskfr(:,3,2,:));         % Lateral central
      fmasklp = squeeze(maskfr(:,1,1,:)&maskfr(:,2,2,:));  % Lateral posterior
%
      fmaskl = {fmasklt; fmasklc; fmasklp};      % Lateral femur
%
      fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));  % Medial trochlea
      fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                       maskfr(:,4,2,:));         % Medial central
      fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));  % Medial posterior
%
      fmaskm = {fmaskmt; fmaskmc; fmaskmp};      % Medial femur
%
% Get Tibia ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
      tmaskla = squeeze(masktr(:,1,1,:));        % Lateral anterior
      tmasklc = squeeze(masktr(:,1,2,:)&masktr(:,2,1,:));  % Lateral central
      tmasklp = squeeze(masktr(:,2,2,:));        % Lateral posterior
%
      tmaskl = {tmaskla; tmasklc; tmasklp};      % Lateral tibia
%
      tmaskma = squeeze(masktr(:,3,1,:));        % Medial anterior
      tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));  % Medial central
      tmaskmp = squeeze(masktr(:,4,2,:));        % Medial posterior
%
      tmaskm = {tmaskma; tmaskmc; tmaskmp};      % Medial tibia
%
% Combine Masks into Cell Arrays
%
      masklay = {maskf; maskt};        % Combine femur and tibia masks
%
      fmask = {fmaskl; fmaskm};        % Combine femur masks
      tmask = {tmaskl; tmaskm};        % Combine tibia masks
      maskroi = {fmask; tmask};        % Combine femur and tibia masks
%
      rslbs = {rslf; rslt};            % Combine femur and tibia slices
%
% Get Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Bone - 1 = femur and 2 = tibia
%   Index 4 - Compartment - 1 = lateral and 2 = medial
%   Index 5 - ROI - 1 = anterior/trochlea, 2 = central and 3 = posterior
%   Index 6 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and compartment analysis variables are:
%        1 = superficial and 2 = deep.
%
      tcp = squeeze(t2s_respx(kss,ileg+1,:,:,:,[2;1]));
      nps = squeeze(t2s_nps(kss,ileg+1,:,:,:,[2;1]));
%
% Plot Results
%
%       sid = [restxt legtxt];      % Subject, result (T1rho/T2*), and leg
      psl_plts(v,masklay,maskroi,rsl,nrsl,rslbs,idt,tcp,nps,mxts, ...
              cmap,sls(ileg+1),pngnamf);
%
   end                  % End of km loop - T2* MAT file loop
%
   close all;           % Close all plot windows
%
end                     % End of ks loop - subjects loop
%
return