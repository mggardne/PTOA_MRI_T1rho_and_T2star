%#######################################################################
%
%                       * DICOM CHecK 3 Program *
%
%          M-File which reads a series of DICOM directories, collects
%     information about the DICOM files, plots the T1 or T2 data, and
%     fits either spin lock (T1rho) or echo times (T2*) for a couple
%     of slices.  Plots are output to Postscript files.
%
%     NOTES:  1.  Reads Philips MRI DICOM files.  May not work with
%             other types of images.
%
%             2.  Program reads and collects information from the first
%             DICOM file header in each series.
%
%             3.  The program traps for some, but not all parameters in
%             the DICOM file header.
%
%             4.  The program looks for subdirectories of the current
%             directory that begin with the letter "s".  The user
%             than selects the subdirectories for processing.
%
%             5.  The program uses the series description to find the
%             spin lock times for the T1rho series.  If the delimiter,
%             "SL", is not found, a warning is issued and the spin lock
%             times are assumed to be 0, 10, 40 and 80 ms.  
%
%             6.  The M-files T1r3d_calc.m, dftreg.m, dftregistration.m
%             and exp1_fun.m must be in the current path or directory.
%
%             7.  The 2D translation image registration is done by the
%             Mathworks File Exchange function, dftregistration.m.  See:
%
% https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation
%
%             8.  The Matlab Optimization toolbox is required.  The
%             nonlinear least squares is performed by the Matlab
%             function lsqcurvefit in the optimization toolbox.
%
%             9.  The nonlinear least squares uses the Levenberg-
%             Marquardt algorithm in the Matlab function lsqcurvefit.
%
%             10.  The Matlab Parallel Computing toolbox is required.
%             The Matlab parallel construct parfor is used to calculate
%             the T1rho/T2* values in parallel.  Use Matlab command
%             parpool to control the number of workers.
%
%     31-Mar-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Series Directories with DICOM Images
%
ddirs = dir('s*');
ddirs = ddirs([ddirs.isdir]');
ddirs = {ddirs.name}';
[isel,ok] = listdlg('ListString',ddirs,'SelectionMode','multiple', ...
                    'Name','MRI Series Directories','PromptString', ...
                    'Select MRI Series: ','ListSize',[150 400]);
%
if ok<1
  error(' *** dicom_chk:  No series selected!');
end
%
ddirs = ddirs(isel);
ns = size(ddirs,1);     % Number of series
%
% Loop through the Series and Get Header Information
%
nimages = zeros(ns,1);  % Number of DICOM files
afiles = cell(ns,1);    % All DICOM file in series
%
etn = zeros(ns,1);      % Echo times as numbers for UTE T2 star sequences
psz = zeros(ns,2);      % Rows and columns
pspc = zeros(ns,2);     % Pixel Spacing
stxt = cell(ns,1);      % Series Description
sn = zeros(ns,1);       % Series number
sdat = datetime(zeros(ns,1),0,0,'Format','dd-MMM-yyyy HH:mm:ss');
%
for k = 1:ns
%
% Get DICOM Files
%
   ddir = ddirs{k};
   dfiles = dir(fullfile(ddir,'i*.dcm'));
   dfiles = dfiles(~[dfiles.isdir]');
   nimages(k) = size(dfiles,1);
   afiles{k} = {dfiles.name}';
%
   fnam = afiles{k}{1};
   fnam = fullfile(ddir,fnam);
%
% Get DICOM Header Information from Files
%
   if exist(fnam,'file')
%
     info = dicominfo(fnam);
%
     if isfield(info,'Rows')
       psz(k,:) = [info.Rows info.Columns];
     end
     if isfield(info,'PixelSpacing')
       pspc(k,:) = info.PixelSpacing';
     end
     stxt{k} = info.SeriesDescription;
     sn(k) = info.SeriesNumber;
     if isfield(info,'EchoTrainLength')
       n = info.EchoTrainLength;
       if n>nimages(k)
         n = nimages(k);
       end
       et = NaN(n,1); % Echo times for this sequence
       et(1) = info.EchoTime;
       for m = 2:n
          fnam = afiles{k}{m};
          fnam = fullfile(ddir,fnam);
          if exist(fnam,'file')
            info1 = dicominfo(fnam);
            et(m) = info1.EchoTime;
          else
            break;
          end
       end
       etu = unique(et);
       idnnan = ~isnan(etu);
       etn(k) = etu(idnnan);
     end
     dat = info.SeriesDate;
     tim = info.SeriesTime;
     sdat(k) = datetime([dat tim],'InputFormat','yyyyMMddHHmmss.SSSSS');
   end
%
end
%
% Get Spin Lock Times from Series Description
%
if contains(stxt,'rho','IgnoreCase',true)
  idvr = contains(stxt,'SL');          % Should contain TSL
%
  if ~any(idvr)
    slt = [0 10 40 80]';
    nslt = 4;
    warning([' *** dicom_chk3:  Not able to get spin lock times', ...
             ' from series description!']);
  else
    slt = extractBetween(stxt(idvr),'SL','ms');  % Spin lock times as text
    slt = strrep(slt,'_',' ');
    slt = strtrim(slt);
    slt = strrep(slt,' ',',');
    slt = eval(['[' slt{1} ']'])';
    nslt = size(slt,1);
    if nslt~=4
      error(' *** dicom_chk3:  Incorrect number of spin lock times!');
    end
  end
end
%
% T1rho
%
if ns==1&&exist('slt','var')           % T1rho
%
  t0 = 65;              % Initial value for T1rho
  mxtr = 80;            % Maximum scale on T1rho plots
%
% Get Index to Files and Get File Names
%
  ddir = ddirs{1};
  fnams = afiles{1};
  nf = size(fnams,1);   % Number of files
  nsl = nf/nslt;        % Number of slices
%
  if sn(1)==1401
    plt_sl = (11:61)';
  else
    plt_sl = floor(nsl/4);
%     plt_sl = plt_sl:plt_sl:3*plt_sl;   % Slices to plot
    plt_sl = [plt_sl; 3*plt_sl];       % Slices to plot
  end
%
  nplts = size(plt_sl,1);
%
%  Get File Name Description and PS File Name
%
  sn1 = int2str(sn(1));
  fs = ['S' sn1];
  ds = [' for Series ' sn1];           % Description with series number
%
  pnam = ['T1rho_' fs '.ps'];          % Plot file name
%
% Set Up Array for Loop
%
  psz = psz(1,:);
  npx = psz(1);
  npy = psz(2);
%
  dat3d = zeros([psz nslt]);
  T1rhonls = zeros(npx,npy,nplts);     % Nonlinear least squares time constant
%
% Loop through Plot Slices
%
  for k = 1:nplts
%
     ks = plt_sl(k);
     sll = int2str(ks);                % Slice number as letters
%
% Loop through Spin Lock Times
%
     for l = 1:nslt     % Loop through spin lock times
%
        n = nslt*(ks-1)+l;
        fnam = fnams{n};          % Filename for this spin lock time
        fprintf(1,['\n Processing file:  ' strrep(fnam,'\','\\') ...
                   ', Slice:  ' sll ', Spin lock time:  ' ...
                   int2str(slt(l)) ' ms']);
%
% Load and Scale Image
%
        img = dicomread(fullfile(ddir,fnam));
        info = dicominfo(fullfile(ddir,fnam));
        sl = double(info.RescaleSlope);
        offst = double(info.RescaleIntercept);   % Usually zero
        r = sl*double(img)+offst;
        dat3d(:,:,l) = r;
%
% Registration of Images for the Different Spin Lock Times
%
        if l>1
          dat3d(:,:,l) = dftreg(dat3d(:,:,l-1),dat3d(:,:,l),100);
        end
     end
%
     fprintf(1,'\n');     % Line between slices
%
% Plot a Montage of Each Slice
%
     figure;
     orient landscape;
     montage(dat3d(:,:,:),'DisplayRange',[0 0.95*max(dat3d(:))]);
     title({['T1' ds]; ['Slice ' sll]},'FontSize',20, ...
             'FontWeight','bold');
     pos = get(gca,'pos');
     pos(2) = 0.02;
     set(gca,'pos',pos);
     if k==1
       print('-dpsc2','-r600','-fillpage',pnam);
     else
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
%
% Calculate T1rho Values for the Whole Slice
%
     T1rhonls(:,:,k) = T1r3d_calc(dat3d,slt,t0);
%
% Plot T1rho Values for this Slice
%
     figure;
     orient landscape;
%
     img = T1rhonls(:,:,k);
     img(img<0) = 0;
     img(img>mxtr) = mxtr;
     imagesc(img,[0 mxtr]);
%
     colormap jet;
     axis image;
     axis off;
%
     colorbar;
%
     title({['T1\rho' ds]; ['Slice ' sll]},'FontSize',16, ...
             'FontWeight','bold');
%
     print('-dpsc2','-r600','-fillpage','-append',pnam);
%
  end
%
% Save T1rho to MAT File
%
  mnam = ['T1rho_' fs '_chk.mat'];     % MAT file name
  save(mnam,'plt_sl','T1rhonls');
%
% T2*
%
elseif ns==5||ns==6                    % T2*
%
  t0 = 35;              % Initial value for T2*
  mxts = 75;            % Maximum scale on T2* plots
%
% Sort Echo Times
%
  [etn,ids] = sort(etn);
  ddirs = ddirs(ids);
  afiles = afiles(ids);
%
% Get Number of Files and Slice Numbers
%
  nf = zeros(ns,1);
%
  for k = 1:ns
     nf(k) = size(afiles{k},1);        % Number of files
  end
%
  if any(diff(nf))
    error([' *** dicom_chk:  Inconsistent number of files in', ...
           ' directories!']);
  end
%
  nsl = nf(1);          % Number of slices = number of files in directories
%
  if sn(1)==501
    plt_sl = (7:41)';
  else
    plt_sl = floor(nsl/4);
%     plt_sl = plt_sl:plt_sl:3*plt_sl;   % Slices to plot
    plt_sl = [plt_sl; 3*plt_sl];       % Slices to plot
  end
%
  nplts = size(plt_sl,1);
%
%  Get File Name Description and PS File Name
%
  sns = int2str(sn);
  sns = [sns repmat(',',ns,1)];
  sna = sns(1,:);
%
  for k = 2:ns
     sna = [sna sns(k,:)];
  end
%
  sna = sna(1:end-1);
  sna = sna(~isspace(sna));
  sna = strrep(sna,',',', ');
%  
  fs = ['S' strip(sns(1,1:end-1)) '_' strip(sns(ns,1:end-1))];
  ds = [' for Series ' sna];           % Description with series number
%
  pnam = ['T2star_' fs '.ps'];         % Plot file name
%
% Set Up Array for Loop
%
  if any(any(diff(psz)))
    error([' *** dicom_chk:  Inconsistent size of images in', ...
           ' directories!']);
  end
%
  psz = psz(1,:);
  npx = psz(1);
  npy = psz(2);
%
  dat3d = zeros([psz ns]);
  T2starnls = zeros(npx,npy,nplts);    % Nonlinear least squares time constant
%
% Loop through Plot Slices
%
  for k = 1:nplts
%   for k = 1:1
%
     ks = plt_sl(k);
     sll = int2str(ks);                % Slice number as letters
%
% Loop through Echo Times/Subdirectories
%
     for l = 1:ns       % Loop through echo times/directories
%
        ddir = ddirs{l};
        fnam = afiles{l}{ks};     % Filename for this echo time and slice
        fnam = fullfile(ddir,fnam);
        fprintf(1,['\n Processing file:  ' strrep(fnam,'\','\\') ...
                   ', Slice:  ' sll ', Echo time:  ' ...
                   num2str(etn(l)) ' ms']);
%
% Load and Scale Image
%
        img = dicomread(fnam);
        info = dicominfo(fnam);
        img = double(img);
        sl = double(info.RescaleSlope);
        offst = double(info.RescaleIntercept);   % Usually zero
        r = sl*double(img)+offst;
        dat3d(:,:,l) = r;
%
% Registration of Images for the Different Echo Times
%
        if l>1
          dat3d(:,:,l) = dftreg(dat3d(:,:,l-1),dat3d(:,:,l),100);
        end
     end
%
     fprintf(1,'\n');     % Line between slices
%
% Plot a Montage of Each Slice
%
     figure;
     orient landscape;
     montage(dat3d(:,:,:),'DisplayRange',[0 0.95*max(dat3d(:))]);
     title({['T2' ds]; ['Slice ' sll]},'FontSize',20, ...
             'FontWeight','bold');
     pos = get(gca,'pos');
     pos(2) = 0.02;
     set(gca,'pos',pos);
     if k==1
       print('-dpsc2','-r600','-fillpage',pnam);
     else
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
%
% Calculate T2* Values for the Whole Slice
%
     T2starnls(:,:,k) = T1r3d_calc(dat3d,etn,t0);
%
% Plot T2* Values for this Slice
%
     figure;
     orient landscape;
%
     img = T2starnls(:,:,k);
     img(img<0) = 0;
     img(img>mxts) = mxts;
     imagesc(img,[0 mxts]);
%
     colormap jet;
     axis image;
     axis off;
%
     colorbar;
%
     title({['T2*' ds]; ['Slice ' sll]},'FontSize',16, ...
             'FontWeight','bold');
%
     print('-dpsc2','-r600','-fillpage','-append',pnam);
%
  end
%
% Save T2* to MAT File
%
  mnam = ['T2star_' fs '_chk.mat'];    % MAT file name
  save(mnam,'dat3d','etn','plt_sl','T2starnls');
%
% Incorrect Number of Directories
%
else
  error(' *** dicom_chk:  Incorrect number of directories!');
end
%
return