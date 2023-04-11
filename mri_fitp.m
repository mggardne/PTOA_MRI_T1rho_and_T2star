%#######################################################################
%
%                       * MRI FIT PTOA Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     MAT files and fits a monoexponential to the MRI data as a function
%     of spin lock or echo times where T1rho or T2* are the time
%     constants of the fits.  Resulting T1rho and T2* values and summary
%     statistics are written to the MS-Excel spreadsheet,
%     mri_fitp.xlsx, in the "Results\mri_fitp" directory.
%
%     NOTES:  1.  Data MAT files must be in subject directories starting
%             with "0".
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "prois".  See rd_m_dicom.m and
%             seg_prois.m.
%
%             3.  M-file exp_fun1.m, cmprt_ana.m and cmprt_plt.m must
%             be in the current directory or path.
%
%     17-Feb-2022 * Mack Gardner-Morse
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
tsmx = 100;             % Maximum valid T2* result
tsmn = 0;               % Minimum valid T2* result
%
mxtr = 80;              % Maximum scale on T1rho plots
mxts = 75;              % Maximum scale on T2* plots
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','mri_fitp');         % Results directory
%
ifirst = true;          % First write to file
xlsnam = 'mri_fitp.xlsx';              % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs1 = {'Subject' 'Result' 'Leg' 'Comprt' 'Bone' ...
         'Layer'};
hdrs2 = {'Pixels' 'T1R/T2S' 'RSS' 'ValidPix' 'Mean' 'Min' 'Max' ...
         'SD' 'COV'};
%
psnam = fullfile(resdir,'mri_fitp_'); % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Subject Directories
%
sdirs = dir('0*');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
nsubj = size(sdirs,1);
%
% Initialize Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Compartment - 1 = lateral and 2 = medial
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 6 - Layer - 1 = deep and 2 = superficial
%
t1r_res = zeros(nsubj,2,2,2,3,2);
t1r_npx = zeros(nsubj,2,2,2,3,2);
t1r_rss = zeros(nsubj,2,2,2,3,2);
%
t1r_respx = cell(nsubj,2,2,2,3,2);
t1r_rsspx = cell(nsubj,2,2,2,3,2);
t1r_nps = cell(nsubj,2,2,2,3,2);
%
t2s_res = zeros(nsubj,2,2,2,3,2);
t2s_npx = zeros(nsubj,2,2,2,3,2);
t2s_rss = zeros(nsubj,2,2,2,3,2);
%
t2s_respx = cell(nsubj,2,2,2,3,2);
t2s_rsspx = cell(nsubj,2,2,2,3,2);
t2s_nps = cell(nsubj,2,2,2,3,2);
%
% Loop through Subjects
%
% for ks = 1:nsubj
for ks = 1:1
%
% Get Subject Directory (Name) and Number
%
   sdir = sdirs{ks};    % Current subject directory (name)
   subj = eval(sdir);   % Subject number
%
   psnams = [psnam sdir];           % Add subject to PS file name
%
% Get T1rho MAT Files in Directory
%
%    ido = false;         % Skip T1rho
   ido = true;          % Do T1rho
%
   if ido
%
     d = dir(fullfile(sdir,'T1rho_S*.mat'));
     roinams = {d.name}';
     idr = contains(roinams,'proi','IgnoreCase',true);     % Masks
%
     rhonams = roinams(~idr);          % Image MAT files
     nrho = size(rhonams,1);
%
     roinams = roinams(idr);           % ROI MAT files
     nroi = size(roinams,1);
%
     if nrho~=nroi
       error([' *** ERROR in mri_fitp:  Number of T1rho MAT files', ...
              ' does not match the number of ROI MAT files!']);
     end
     clear nroi;
%
% T1rho Identifier
%
     ires = 0;          % ires = 0 - T1rho, ires = 1 - T2*
     idt = 1;           % Spin lock/echo time for plots - 1 = 0 ms spin lock time
%
     psnamr = [psnams '_T1R_'];        % Add result type to PS file name
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
          ileg = 0;     % Coding for leg
        else
          leg = 'R';
          ileg = 1;
        end
%
% Add Leg to PS File Name
%
        psnamf = [psnamr leg '_' pstyp];    % Add leg and load to PS file name
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
        idfl = cell2mat(cellfun(@any,fmaskl,'UniformOutput',false));
        idfl = any(idfl);
        rslfl = rslf(idfl);
%
        fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));     % Medial trochlea
        fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                         maskfr(:,4,2,:));       % Medial central
        fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));     % Medial posterior
%
        fmaskm = {fmaskmt; fmaskmc; fmaskmp};    % Medial femur
        idfm = cell2mat(cellfun(@any,fmaskm,'UniformOutput',false));
        idfm = any(idfm);
        rslfm = rslf(idfm);
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
        idtl = cell2mat(cellfun(@any,tmaskl,'UniformOutput',false));
        idtl = any(idtl);
        rsltl = rslt(idtl);
%
        tmaskma = squeeze(masktr(:,3,1,:));      % Medial anterior
        tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));     % Medial central
        tmaskmp = squeeze(masktr(:,4,2,:));      % Medial posterior
%
        tmaskm = {tmaskma; tmaskmc; tmaskmp};    % Medial tibia
        idtm = cell2mat(cellfun(@any,tmaskm,'UniformOutput',false));
        idtm = any(idtm);
        rsltm = rslt(idtm);
%
% Combine Masks into Cell Arrays
%
        masklay = {maskf; maskt};      % Combine femur and tibia masks
%
        maskl = {fmaskl; tmaskl};      % Combine femur and tibia masks
        maskm = {fmaskm; tmaskm};      % Combine femur and tibia masks
        maskroi = {maskl; maskm};      % Combine compartment masks
%
        rslbs = {rslf; rslt};          % Combine femur and tibia slices
%
        rsll = union(rslfl,rsltl);
        rslm = union(rslfm,rsltm);
        rsls = {rsll; rslm};           % Lateral - row 1, medial - row 2
        nrsls = [size(rsll,1); size(rslm,1)];
%
% Do Compartmental Analysis
%
keyboard
        [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = pcmprt_ana(v,masklay, ...
               maskroi,rsls,nrsls,rsl,rslbs,splt,nslt,fun,init,tr0,opt);
        na = size(tc,1);               % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Compartment - 1 = lateral and 2 = medial
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 6 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and compartment analysis variables are:
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
        sid = ['Subject ' sdir];
        pcmprt_plt(v,masklay,maskroi,rsls,nrsls,rsl,rslbs,idt,tcp, ...
                   nps,mxtr,cmap,sid,psnamf);
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
           tcpv = tcp{ka}(idv);        % Valid T1rho values
           tcpm(ka) = mean(tcpv);      % Mean
           tcpmn(ka) = min(tcpv);      % Minimum
           tcpmx(ka) = max(tcpv);      % Maximum
           tcpsd(ka) = std(tcpv);      % SD
        end
%
        tcpcov = 100*tcpsd./tcpm;      % Coefficient of variation
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
   starnams = roinams(~idr);           % Image MAT files
   nstar = size(starnams,1);
%
   roinams = roinams(idr);             % ROI MAT files
   nroi = size(roinams,1);
%
   if nstar~=nroi
     error([' *** ERROR in mri_fitp:  Number of T2* MAT files', ...
            ' does not match the number of ROI MAT files!']);
   end
   clear nroi;
%
% T2* Identifier
%
   ires = 1;            % ires = 0 - T1rho, ires = 1 - T2*
   idt = 3;             % Spin lock/echo time for plots - 3 = 5 ms echo time
%
   psnamr = [psnamv '_T2S_'];          % Add result type to PS file name
%
% Loop through T2* MAT Files
%
   for km = 1:nstar
%
% Load Data
%
      starnam = starnams{km};
      load(fullfile(sdir,starnam),'etns','iszs','netn','scmx', ...
           'sns','snt','st','v');
      npix = prod(iszs);     % Number of pixels in an image
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
        ileg = 0;       % Coding for leg
      else
        leg = 'R';
        ileg = 1;
      end
%
% Add Leg to PS File Name
%
      psnamf = [psnamr leg '_' pstyp];      % Add leg and load to PS file name
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
      idfl = cell2mat(cellfun(@any,fmaskl,'UniformOutput',false));
      idfl = any(idfl);
      rslfl = rslf(idfl);
%
      fmaskmt = squeeze(maskfr(:,1,2,:)&maskfr(:,4,1,:));  % Medial trochlea
      fmaskmc = squeeze(maskfr(:,1,2,:)&maskfr(:,2,1,:)& ...
                       maskfr(:,4,2,:));         % Medial central
      fmaskmp = squeeze(maskfr(:,1,2,:)&maskfr(:,2,2,:));  % Medial posterior
%
      fmaskm = {fmaskmt; fmaskmc; fmaskmp};      % Medial femur
      idfm = cell2mat(cellfun(@any,fmaskm,'UniformOutput',false));
      idfm = any(idfm);
      rslfm = rslf(idfm);
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
      idtl = cell2mat(cellfun(@any,tmaskl,'UniformOutput',false));
      idtl = any(idtl);
      rsltl = rslt(idtl);
%
      tmaskma = squeeze(masktr(:,3,1,:));        % Medial anterior
      tmaskmc = squeeze(masktr(:,3,2,:)&masktr(:,4,1,:));  % Medial central
      tmaskmp = squeeze(masktr(:,4,2,:));        % Medial posterior
%
      tmaskm = {tmaskma; tmaskmc; tmaskmp};      % Medial tibia
      idtm = cell2mat(cellfun(@any,tmaskm,'UniformOutput',false));
      idtm = any(idtm);
      rsltm = rslt(idtm);
%
% Combine Masks into Cell Arrays
%
      masklay = {maskf; maskt};        % Combine femur and tibia masks
%
      maskl = {fmaskl; tmaskl};        % Combine femur and tibia masks
      maskm = {fmaskm; tmaskm};        % Combine femur and tibia masks
      maskroi = {maskl; maskm};        % Combine compartment masks
%
      rslbs = {rslf; rslt};            % Combine femur and tibia slices
%
      rsll = union(rslfl,rsltl);
      rslm = union(rslfm,rsltm);
      rsls = {rsll; rslm};             % Lateral - row 1, medial - row 2
      nrsls = [size(rsll,1); size(rslm,1)];
%
% Do Compartmental Analysis
%
      [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = pcmprt_ana(v,masklay, ...
               maskroi,rsls,nrsls,rsl,rslbs,splt,nslt,fun,init,tr0,opt);
      na = size(tc,1);              % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Compartment - 1 = lateral and 2 = medial
%   Index 4 - Bone - 1 = femur and 2 = tibia
%   Index 5 - ROI - 1 = anterior/trochlea, 2 = central and 3 - posterior
%   Index 6 - Layer - 1 = deep and 2 = superficial
%
% Note:  Layers for masks and compartment analysis variables are:
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
      sid = ['Subject ' sdir];
      pcmprt_plt(v,masklay,maskroi,rsls,nrsls,rsl,rslbs,idt,tcp, ...
                 nps,mxtr,cmap,sid,psnamf);
%
% Get Statistics on Pixel Results
%
      npxv = zeros(na,1);              % Number of valid results
      tcpm = zeros(na,1);              % Mean
      tcpmn = zeros(na,1);             % Minimum
      tcpmx = zeros(na,1);             % Maximum
      tcpsd = zeros(na,1);             % SD
%
      for ka = 1:na
         idv = tcp{ka}>=tsmn&tcp{ka}<=tsmx;
         npxv(ka) = sum(idv);          % Number of valid results
         tcpv = tcp{ka}(idv);          % Valid T2* values
         tcpm(ka) = mean(tcpv);        % Mean
         tcpmn(ka) = min(tcpv);        % Minimum
         tcpmx(ka) = max(tcpv);        % Maximum
         tcpsd(ka) = std(tcpv);        % SD
      end
%
      tcpcov = 100*tcpsd./tcpm;        % Coefficient of variation
%
% Combine Identifiers
%
      ids = [subj ires ileg];          % MAT file identifiers
      ids = repmat(ids,na,1);
      ids = [ids id];                  % All identifiers
%
% Create and Write Table of Results
%
      t1 = array2table(ids,'VariableNames',hdrs1);
      t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                 'VariableNames',hdrs2);
      t = [t1 t2];
%
      writetable(t,xlsnam,'WriteMode','append', ...
                 'WriteVariableNames',false);
%
   end                  % End of km loop - T2* MAT file loop
%
   close all;           % Close all plot windows
%
end                     % End of ks loop - subjects loop
%
% Save to MAT File
%
save(fullfile(resdir,'mri_fitp.mat'),'t1r_res','t1r_npx','t1r_rss', ...
     't1r_respx','t1r_rsspx','t1r_nps','t2s_res','t2s_npx', ...
     't2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
%
return