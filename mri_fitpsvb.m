%#######################################################################
%
%             * MRI FIT PTOA Slice Visit Baseline Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     MAT files and fits a monoexponential to the MRI data as a function
%     of spin lock or echo times where T1rho or T2* are the time
%     constants of the fits.  Resulting T1rho and T2* values and summary
%     statistics are written to the MS-Excel spreadsheet,
%     mri_fitpsvb.xlsx, in the "Results\mri_fitpsvb" directory.
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "0".
%
%             2.  For BASELINE visit only.
%
%             3.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "prois".  See rd_m_dicom.m and
%             seg_prois_b.m.
%
%             4.  M-file exp_fun1.m, psl_ana.m and psl_plt.m must be in
%             the current directory or path.
%
%     30-Apr-2025 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Map
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter', ...
               2e+3,'Algorithm','levenberg-marquardt','Jacobian', ...
               'on','UseParallel',true);
%
fun = @exp_fun1;        % Exponential function
%
% Initialize Parameters
%
% init = -1;              % Use weighted least squares for starting parameters
% init = 0;               % Use linear least squares for starting parameters
init = 1;               % Use fixed starting parameters
tr0 = 65;               % Initial T1rho estimate in ms
% tr0 = 80;               % Initial T1rho estimate in ms
trmx = 100;             % Maximum valid T1rho result
trmn = 0;               % Minimum valid T1rho result
ts0 = 35;               % Initial T2* estimate in ms
tsmx = 65;              % Maximum valid T2* result
tsmn = 0;               % Minimum valid T2* result
%
mxtr = 80;              % Maximum scale on T1rho plots
mxts = 65;              % Maximum scale on T2* plots
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','mri_fitpsvb');      % Results directory
%
ifirst = true;          % First write to file
% ifirst = false;         % First write to file
xlsnam = 'mri_fitpsvb.xlsx';           % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs = {'Subject' 'Visit' 'Result' 'Leg' 'Bone' 'Comprt' 'ROI' 'Layer'};
hdrs2 = {'Pixels' 'T1R/T2S' 'RSS' 'ValidPix' 'Mean' 'Min' 'Max' ...
         'SD' 'COV'};
%
psnam = fullfile(resdir,'mri_fitpsvb_');         % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Subject and Visit Directories
%
sdirs = dir('0*');
sdirs = {sdirs([sdirs.isdir]').name}';         % Subject directories
nsubj = size(sdirs,1);
% subjn = [(5:8)';15;17;20;21;23;27];
% nsubj = size(subjn,1);
%
% sdirs = int2str(subjn);
% sdirs = [repmat('0',nsubj,1) sdirs];
% sdirs = strrep(cellstr(sdirs),' ','0');
vdirs = {'BASELINE';'YEAR1';'YEAR2'};
nvdirs = size(vdirs,1);
vtxts = ['BL';'Y1';'Y2'];              % Abbreviations of visits
%
% Initialize Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit - 1 = BASELINE, 2 = YEAR1, and 3 = YEAR2
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 7 - Layer - 1 = deep and 2 = superficial
%
t1r_res = zeros(nsubj,3,2,2,2,3,2);
t1r_npx = zeros(nsubj,3,2,2,2,3,2);
t1r_rss = zeros(nsubj,3,2,2,2,3,2);
%
t1r_respx = cell(nsubj,3,2,2,2,3,2);
t1r_rsspx = cell(nsubj,3,2,2,2,3,2);
t1r_nps = cell(nsubj,3,2,2,2,3,2);
%
t2s_res = zeros(nsubj,3,2,2,2,3,2);
t2s_npx = zeros(nsubj,3,2,2,2,3,2);
t2s_rss = zeros(nsubj,3,2,2,2,3,2);
%
t2s_respx = cell(nsubj,3,2,2,2,3,2);
t2s_rsspx = cell(nsubj,3,2,2,2,3,2);
t2s_nps = cell(nsubj,3,2,2,2,3,2);
%
% If NOT First Run, Load Previous Results
%
if ~ifirst
  load(fullfile(resdir,'mri_fitpsvb.mat'));
end
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory (Name) and Number
%
   sdir = sdirs{ks};    % Current subject directory (name)
%    subj = subjn(ks);    % Subject number
   subj = str2num(sdir);               % Subject number
   subjtxt = ['Subject ' sdir];
%
   psnams = [psnam sdir];           % Add subject to PS file name
%
% Loop through Visits
%
%    for kv = 1:nvdirs
   for kv = 1           % Baseline only
%
      vdir = vdirs{kv}; % Current visit directory
      vtxt = vtxts(kv,:);
%
      psnamsv = [psnams '_' vtxt];     % Add visit to PS file name
%
% Visit Identifier
%
      ivis = kv-1;      % ivis - 0 = BASELINE, 1 = YEAR1, and 2 = YEAR2
%
% Get T1rho MAT Files in Directory
%
%       ido = false;      % Skip T1rho
      ido = true;       % Do T1rho
%
      if ido
%
        d = dir(fullfile(sdir,vdir,'T1rho_S*.mat'));
        roinams = {d.name}';
        idr = contains(roinams,'01_prois','IgnoreCase',true);   % Masks
%
        rhonams = roinams(~idr);
        idx = contains(rhonams,'chk','IgnoreCase',true);   % Check files
        rhonams = rhonams(~idx);
        idx = contains(rhonams,'mrois','IgnoreCase',true); % Meniscus ROI files
        rhonams = rhonams(~idx);
        idx = contains(rhonams,'2dprois','IgnoreCase',true);    % 2D ROI files
        rhonams = rhonams(~idx);       % Image MAT files
        nrho = size(rhonams,1);
%
        roinams = roinams(idr);        % 3D ROI MAT files
        nroi = size(roinams,1);
%
        if nrho~=nroi
          error([' *** ERROR in mri_fitpsvb:  Number of T1rho MAT', ...
                 ' files does not match the number of ROI MAT files!']);
        end
        clear nroi;
%
% T1rho Identifier
%
        ires = 0;       % ires = 0 - T1rho, ires = 1 - T2*
        restxt = [subjtxt ' ' vtxt ' T1\rho'];
        idt = 1;        % Spin lock/echo time for plots - 1 = 0 ms spin lock time
%
        psnamr = [psnamsv '_T1R_'];    % Add result type to PS file name
%
% Loop through T1rho MAT Files
%
        for km = 1:nrho
%
% Load Data
%
           rhonam = rhonams{km};
           load(fullfile(sdir,vdir,rhonam),'iszs','nslt','scmx', ...
                'sns','snt','splt','st','v');
           npix = prod(iszs);          % Number of pixels in an image
           fs = ['S' snt];        % Series number prefaced with a 'S'
%
           idm = contains(roinams,rhonam(1:end-4));   % Get matching file
           roinam = roinams{idm};
           load(fullfile(sdir,vdir,roinam),'maskf','maskfr', ...
                'maskt','masktr','nrsl','rsl','rslf','rslt');
%
% Parse Series Text for Leg
%
           if strcmpi(st(1),'L')
             leg = 'L';
             legtxt = ' Left Leg';
             ileg = 0;  % Coding for leg
           else
             leg = 'R';
             legtxt = ' Right Leg';
             ileg = 1;  % Coding for leg
           end
%
% Add Leg to PS File Name
%
           psnamf = [psnamr leg pstyp];     % Add leg to PS file name
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
                            maskfr(:,4,2,:));       % Medial central
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
           tmaskla = squeeze(masktr(:,1,1,:));   % Lateral anterior
           tmasklc = squeeze(masktr(:,1,2,:)&masktr(:,2,1,:));  % Lateral central
           tmasklp = squeeze(masktr(:,2,2,:));   % Lateral posterior
%
           tmaskl = {tmaskla; tmasklc; tmasklp};      % Lateral tibia
%
           tmaskma = squeeze(masktr(:,3,1,:));   % Medial anterior
           tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));  % Medial central
           tmaskmp = squeeze(masktr(:,4,2,:));   % Medial posterior
%
           tmaskm = {tmaskma; tmaskmc; tmaskmp};      % Medial tibia
%
% Combine Masks into Cell Arrays
%
           masklay = {maskf; maskt};   % Combine femur and tibia masks
%
           fmask = {fmaskl; fmaskm};   % Combine femur masks
           tmask = {tmaskl; tmaskm};   % Combine tibia masks
           maskroi = {fmask; tmask};   % Combine femur and tibia masks
%
           rslbs = {rslf; rslt};       % Combine femur and tibia slices
%
% Do Slice Analysis
%
           [tc,~,rss,npx,id,tcp,~,rssp,nps] = psl_ana(v,masklay, ...
                     maskroi,rsl,nrsl,rslbs,splt,nslt,fun,init,tr0,opt);
           na = size(tc,1);            % Number of results
%
% Set Missing Data to NaNs
%
           idnan = npx<1;
           if any(idnan)
             tc(idnan) = NaN;
             tcp(idnan) = {NaN};
           end
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit - 1 = BASELINE, 2 = YEAR1, and 3 = YEAR2
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 7 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and analysis variables are:
%        1 = superficial and 2 = deep.
%
           for ka = 1:na
              t1r_res(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                      id(ka,4)+1) = tc(ka);
              t1r_npx(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                      id(ka,4)+1) = npx(ka);
              t1r_rss(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                      id(ka,4)+1) = rss(ka);
              t1r_respx{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1, ...
                        id(ka,3)+1,id(ka,4)+1} = tcp{ka};
              t1r_rsspx{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1, ...
                        id(ka,3)+1,id(ka,4)+1} = rssp{ka};
              t1r_nps{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                        id(ka,4)+1} = nps{ka};
           end
%
% Plot Results
%
           sid = [restxt legtxt];      % Subject, visit, result (T1rho/T2*), and leg
           psl_plt(v,masklay,maskroi,rsl,nrsl,rslbs,idt,tcp,nps, ...
                   mxtr,cmap,sid,psnamf);
%
% Get Statistics on Pixel Results
%
           npxv = zeros(na,1);         % Number of valid results
           tcpm = zeros(na,1);         % Mean
           tcpmn = zeros(na,1);        % Minimum
           tcpmx = zeros(na,1);        % Maximum
           tcpsd = zeros(na,1);        % SD
%
           for ka = 1:na
              idv = tcp{ka}>=trmn&tcp{ka}<=trmx;
              npxv(ka) = sum(idv);     % Number of valid results
              if npxv(ka)>0
                tcpv = tcp{ka}(idv);   % Valid T1rho values
                tcpm(ka) = mean(tcpv); % Mean
                tcpmn(ka) = min(tcpv); % Minimum
                tcpmx(ka) = max(tcpv); % Maximum
                tcpsd(ka) = std(tcpv); % SD
              end
           end
%
           tcpcov = 100*tcpsd./tcpm;   % Coefficient of variation
           tcpcov(isnan(tcpcov)) = 0;  % Catch any NaNs
%
% Combine Identifiers
%
           ids = [subj ivis ires ileg];     % MAT file identifiers
           ids = repmat(ids,na,1);
           ids = [ids id];             % All identifiers
%
% Create and Write Table of Results
%
           t1 = array2table(ids,'VariableNames',hdrs);
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
        end             % End of km loop - T1rho MAT file loop
%
        close all;      % Close all plot windows
%
      end               % End of ido - Skip T1rho?
%
% Get T2* MAT Files in Directory
%
      d = dir(fullfile(sdir,vdir,'T2star_S*.mat'));
      roinams = {d.name}';
      idr = contains(roinams,'01_prois','IgnoreCase',true);     % Masks
%
      starnams = roinams(~idr);
      idx = contains(starnams,'chk','IgnoreCase',true);    % Check files
      starnams = starnams(~idx);
      idx = contains(starnams,'mrois','IgnoreCase',true);  % Meniscus ROI files
      starnams = starnams(~idx);
      idx = contains(starnams,'2dprois','IgnoreCase',true);     % 2D ROI files
      starnams = starnams(~idx);       % Image MAT files
      nstar = size(starnams,1);
%
      roinams = roinams(idr);          % 3D ROI MAT files
      nroi = size(roinams,1);
%
      if nstar~=nroi
        error([' *** ERROR in mri_fitpsvb:  Number of T2* MAT files', ...
               ' does not match the number of ROI MAT files!']);
      end
      clear nroi;
%
% T2* Identifier
%
      ires = 1;         % ires = 0 - T1rho, ires = 1 - T2*
      restxt = [subjtxt ' ' vtxt ' T2*'];
      idt = 4;          % Spin lock/echo time for plots - 4 = 5 ms echo time
%
      psnamr = [psnamsv '_T2S_'];      % Add result type to PS file name
%
% Loop through T2* MAT Files
%
      for km = 1:nstar
%
% Load Data
%
         starnam = starnams{km};
         load(fullfile(sdir,vdir,starnam),'etns','iszs','netn', ...
              'scmx','sns','snt','st','v');
         npix = prod(iszs);            % Number of pixels in an image
         fs = ['S' snt];          % Series number prefaced with a 'S'
%
         idm = contains(roinams,starnam(1:end-4));    % Get matching file
         roinam = roinams{idm};
         load(fullfile(sdir,vdir,roinam),'maskf','maskfr', ...
              'maskt','masktr','nrsl','rsl','rslf','rslt');
%
% Parse Series Text for Leg
%
         if strcmpi(st(1),'L')
           leg = 'L';
           legtxt = ' Left Leg';
           ileg = 0;    % Coding for leg
         else
           leg = 'R';
           legtxt = ' Right Leg';
           ileg = 1;    % Coding for leg
         end
%
% Add Leg to PS File Name
%
         psnamf = [psnamr leg pstyp];  % Add leg to PS file name
%
% Get Femur ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral-medial, posterior-center, lateral trochlea, and medial trochlea)
%   Third:   Sides of plane (lateral/medial, center/posterior, lateral trochlea/center, and medial trochlea/center)
%   Fourth:  Number of slices
%
         fmasklt = squeeze(maskfr(:,1,1,:)&maskfr(:,3,1,:));    % Lateral trochlea
         fmasklc = squeeze(maskfr(:,1,1,:)&maskfr(:,2,1,:)& ...
                          maskfr(:,3,2,:));      % Lateral central
         fmasklp = squeeze(maskfr(:,1,1,:)&maskfr(:,2,2,:));    % Lateral posterior
%
         fmaskl = {fmasklt; fmasklc; fmasklp};   % Lateral femur
%
         fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));    % Medial trochlea
         fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                          maskfr(:,4,2,:));      % Medial central
         fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));    % Medial posterior
%
         fmaskm = {fmaskmt; fmaskmc; fmaskmp};   % Medial femur
%
% Get Tibia ROI Masks
%
% Dimesions:
%   First:   Number of pixel in slice image
%   Second:  Number of planes (lateral central-anterior, lateral posterior-central, medial central-anterior, medial posterior-central)
%   Third:   Sides of plane (anterior/central, or central/posterior)
%   Fourth:  Number of slices
%
         tmaskla = squeeze(masktr(:,1,1,:));     % Lateral anterior
         tmasklc = squeeze(masktr(:,1,2,:)&masktr(:,2,1,:));    % Lateral central
         tmasklp = squeeze(masktr(:,2,2,:));     % Lateral posterior
%
         tmaskl = {tmaskla; tmasklc; tmasklp};   % Lateral tibia
%
         tmaskma = squeeze(masktr(:,3,1,:));     % Medial anterior
         tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));    % Medial central
         tmaskmp = squeeze(masktr(:,4,2,:));     % Medial posterior
%
         tmaskm = {tmaskma; tmaskmc; tmaskmp};   % Medial tibia
%
% Combine Masks into Cell Arrays
%
         masklay = {maskf; maskt};     % Combine femur and tibia masks
%
         fmask = {fmaskl; fmaskm};     % Combine femur masks
         tmask = {tmaskl; tmaskm};     % Combine tibia masks
         maskroi = {fmask; tmask};     % Combine femur and tibia masks
%
         rslbs = {rslf; rslt};         % Combine femur and tibia slices
%
% Do Slice Analysis
%
         [tc,~,rss,npx,id,tcp,~,rssp,nps] = psl_ana(v,masklay, ...
                        maskroi,rsl,nrsl,rslbs,etns,netn,fun,init,ts0,opt);
         na = size(tc,1);              % Number of results
%
% Set Missing Data to NaNs
%
         idnan = npx<1;
         if any(idnan)
           tc(idnan) = NaN;
           tcp(idnan) = {NaN};
         end
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Visit - 1 = BASELINE, 2 = YEAR1, and 3 = YEAR2
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - Compartment - 1 = lateral and 2 = medial
%   Index 6 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 7 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and compartment analysis variables are:
%        1 = superficial and 2 = deep.
%
         for ka = 1:na
            t2s_res(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                    id(ka,4)+1) = tc(ka);
            t2s_npx(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                    id(ka,4)+1) = npx(ka);
            t2s_rss(ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                    id(ka,4)+1) = rss(ka);
            t2s_respx{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                      id(ka,4)+1} = tcp{ka};
            t2s_rsspx{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                      id(ka,4)+1} = rssp{ka};
            t2s_nps{ks,kv,ileg+1,id(ka,1)+1,id(ka,2)+1,id(ka,3)+1, ...
                      id(ka,4)+1} = nps{ka};
         end
%
% Plot Results
%
         sid = [restxt legtxt];   % Subject, visit, result (T1rho/T2*), and leg
         psl_plt(v,masklay,maskroi,rsl,nrsl,rslbs,idt,tcp,nps,mxts, ...
                 cmap,sid,psnamf);
%
% Get Statistics on Pixel Results
%
         npxv = zeros(na,1);           % Number of valid results
         tcpm = zeros(na,1);           % Mean
         tcpmn = zeros(na,1);          % Minimum
         tcpmx = zeros(na,1);          % Maximum
         tcpsd = zeros(na,1);          % SD
%
         for ka = 1:na
            idv = tcp{ka}>=tsmn&tcp{ka}<=tsmx;
            npxv(ka) = sum(idv);       % Number of valid results
            if npxv(ka)>0
              tcpv = tcp{ka}(idv);     % Valid T2* values
              tcpm(ka) = mean(tcpv);   % Mean
              tcpmn(ka) = min(tcpv);   % Minimum
              tcpmx(ka) = max(tcpv);   % Maximum
              tcpsd(ka) = std(tcpv);   % SD
            end
         end
%
         tcpcov = 100*tcpsd./tcpm;     % Coefficient of variation
         tcpcov(isnan(tcpcov)) = 0;    % Catch any NaNs
%
% Combine Identifiers
%
         ids = [subj ivis ires ileg];  % MAT file identifiers
         ids = repmat(ids,na,1);
         ids = [ids id];               % All identifiers
%
% Create and Write Table of Results
%
         t1 = array2table(ids,'VariableNames',hdrs);
         t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                    'VariableNames',hdrs2);
         t = [t1 t2];
%
         writetable(t,xlsnam,'WriteMode','append', ...
                    'WriteVariableNames',false);
%
      end               % End of km loop - T2* MAT file loop
%
      close all;        % Close all plot windows
%
   end                  % End of kv loop - visits loop
%
end                     % End of ks loop - subjects loop
%
% Save to MAT File
%
save(fullfile(resdir,'mri_fitpsvb.mat'),'t1r_res','t1r_npx', ...
     't1r_rss','t1r_respx','t1r_rsspx','t1r_nps','t2s_res', ...
     't2s_npx','t2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
%
return