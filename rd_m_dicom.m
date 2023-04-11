%#######################################################################
%
%              * ReaD Meniscus DICOM Volume Data Program *
%
%          M-File which reads the meniscus DICOM images from the T1rho
%     and T2* series.  The image data is put into volume matrices and
%     the volumes from different spin lock or echo times are registered
%     to the spin lock time = 0 ms volume for T1rho or echo time = 5 ms
%     volume for T2*.  The registered volume data and series
%     information are saved to separate MAT files.
%
%     NOTES:  1.  Matlab MAT file dicom_lst2.mat must be in the current
%             directory or path.
%
%             2.  elastix.exe must be installed and executable on this
%             computer.
%
%             3.  The Matlab path must include the following two paths:
% C:\Users\mggardne\BRUCE\Risk_Fac\CACL\MRI_data\MelastiX\yamlmatlab-master
% C:\Users\mggardne\BRUCE\Risk_Fac\CACL\MRI_data\MelastiX\matlab_elastix-master\code
%             These are required to run elastix using Rob Campbell's
%             MelastiX.
%
%             4.  Elastix parameter file, Parameters_RigidBody.txt,
%             must in the current directory.
%
%             5.  The M-files dftreg.m and dftregistration.m must be in
%             the current path or directory.
%
%             6.  The Matlab Image Processing Toolbox is required.
%
%             7.  The spin lock and echo times must be certain values
%             and in ascending order.  The last duplicate echo times
%             are used for registration and saved in the resulting MAT
%             file.

%     19-Jul-2022 * Mack Gardner-Morse
%

%#######################################################################
%
rad2deg = 180/pi;       % Radians to degrees
%
% Load Image File Data from dicom_lst2.mat
%
load dicom_lst2.mat afiles ddirs etn ets idvr isz nimages pspc sn ...
                    splcktc stxt;
%
% Use Series Descriptions to Find T1rho Series
%
idr = contains(stxt,'rho','IgnoreCase',true)&nimages>99;   % T1rho
nr = sum(idr);
%
ddirr = ddirs(idr);     % Subdirectories for T1rho series
nfiler = nimages(idr);  % Numbers of T1rho files
afiler = afiles(idr);   % T1rho files
iszr = isz(idr,:);      % Image sizes in pixels
snr = sn(idr);          % Series numbers
pspcr = pspc(idr,:);    % Pixel sizes
spltr = splcktc(idr);   % T1rho spin lock times
stxtr = stxt(idr);      % T1rho series
%
% Check Spin Lock Times
%
ichk = ~contains(spltr,'0,10,40,80')';
%
for k = find(ichk)
   slt = extractBetween(stxtr(k),'0 ',textBoundary,'Boundaries', ...
                       'inclusive');
   slt = strrep(slt,'ms','');
   slt = strrep(slt,'_',' ');
   slt = strtrim(slt);
   slt = strrep(slt,' ',',');
   spltr(k) = slt;
end
%
% Get Spin Lock Times as a Matrix
%
sltm = cell(1,nr);
%
for k = 1:nr
   sltm{k} = eval(['[' spltr{k} ']'])';
   nslt = size(sltm{k},1);
   if nslt~=4
     error(' *** rd_m_dicom:  Incorrect number of spin lock times!');
   end
end
%
sltm = cell2mat(sltm);  % Spin lock times in a matrix
%
% Loop through the T1rho Series
%
nreg = nr*(nslt-1);     % Number of volume registrations
tmr = zeros(nreg,10);   % 3D and 2D registration translations/rotations
sltt = cell(nreg,1);    % Registration spin lock times
sr = cell(nreg,1);      % Series
fitr = string(repmat('T1rho',nreg,1)); % Type of fit (T1rho or T2star)
%
for k = 1:nr
%
% Get Series Specific Indices
%
   ddir = ddirr{k};     % Subdirectory for series
%
   nfile = nfiler(k);   % Number of image files in this series
   splt = eval(['[' spltr{k} ']'])';   % T1rho spin lock times
   nsls = nfile/nslt;   % Number of slices
   iszs = iszr(k,:);    % Image size
   pspcs = pspcr(k,:);  % Pixel size
%
   sns = snr(k);        % Series number
   snt = int2str(sns);  % Series number as text
   psfile = ['S' snt '.ps'];           % PS file name
%
% Get Image File Names
%
   fnams = afiler{k};   % T1rho files
%
% Get Images and Maximum Scaled Image Values
%
   valmx = zeros(nfile,1);
   v = zeros([iszs nsls nslt]);
%
   for l = 1:nsls
      for m = 1:nslt
         idx = nslt*(l-1)+m;
         info = dicominfo(fullfile(ddir,fnams{idx}));
         sl = double(info.RescaleSlope);
         y0 = double(info.RescaleIntercept);
         img = dicomread(info);
         img = sl*double(img)+y0;
         v(:,:,l,m) = img;
         valmx(idx) = max(img(:));     % Slice maximum
      end
   end
%
   scmx = 10*fix(max(valmx)/10);       % Round maximum value down
%
% Register the Different Spin Lock Time Images to the 0 ms Spin Lock
% Time Image
%
   rsl(2,1) = floor(nsls/4);
   rsl = [rsl(2,1); 3*rsl(2,1)];       % Slices for 2D registration
%
   pxmx = 3*iszs(1);                   % Maximum X pixels
   px = 0:32:pxmx;                     % Grid of 32 pixels
   nx = size(px,2);                    % Number of X grid lines
   pymx = iszs(2);                     % Maximum Y pixels
   py = 0:32:pymx;                     % Grid of 32 pixels
   ny = size(py,2);                    % Number of Y grid lines
%
   t2 = cell(2,1);      % 2D translations
%
   for l = 2:nslt
%
      o = nslt*k-nslt-k+l;
      sr{o} = ['S' snt];
      sltt{o} = ['Spin lock time ' int2str(splt(l)),' to ', ...
                  int2str(splt(l-1)) ' ms'];
%
      [reg,t] = elastix(v(:,:,:,l),v(:,:,:,l-1),[], ...
                               'Parameters_RigidBody.txt');     % 3D
%
      t = t.TransformParameters{1};

      rx = t.TransformParameters(1);
      ry = t.TransformParameters(2);
      rz = t.TransformParameters(3);
      tx = t.TransformParameters(4);
      ty = t.TransformParameters(5);
      tz = t.TransformParameters(6);
      tn = [tx; ty; tz; rad2deg*[rx; ry; rz]];
      t = sprintf(['tx = %.1f, ty = %.1f, tz = %.1f, rx =  %.1f, ', ...
                   'ry =  %.1f, rz =  %.1f'],tn);
%
      [tf1,dr1,dc1] = dftreg(v(:,:,rsl(1),l-1), ...
                             v(:,:,rsl(1),l),100);         % 2D
      t2{1} = sprintf('tx = %.1f, ty = %.1f',dc1,dr1);
      [tf2,dr2,dc2] = dftreg(v(:,:,rsl(2),l-1), ...
                             v(:,:,rsl(2),l),100);         % 2D
      t2{2} = sprintf('tx = %.1f, ty = %.1f',dc2,dr2);
%
      tmr(o,:) = [tn',dc1,dr1,dc2,dr2];
%
% Plot Registration
%
      for m = 1:2
%
% Plots with Grids
%
        vp = squeeze(v(:,:,rsl(m),l-1:l));
        irng = [min(vp(:)) max(vp(:))];
        figure;
        subplot(2,1,1);
        montage({vp(:,:,1),vp(:,:,2),reg(:,:,rsl(m))}, ...
                'Size',[1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        orient landscape;
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t,'FontSize',12,'FontWeight','bold');
        ylabel('3D','FontSize',12,'FontWeight','bold');
        title({['Series ' snt]; sltt{o}; ['Slice ', ...
              int2str(rsl(m))]},'FontSize',16,'FontWeight','bold');
%
        subplot(2,1,2);
        if m==1
          montage({vp(:,:,1),vp(:,:,2),tf1},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        else
          montage({vp(:,:,1),vp(:,:,2),tf2},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        end
        colormap gray;
        brighten(l*l/32);
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t2{m},'FontSize',12,'FontWeight','bold');
        ylabel('2D','FontSize',12,'FontWeight','bold');
%
        if l==2&&m==1
          print('-dpsc2','-r600','-fillpage',psfile);
        else
          print('-dpsc2','-r600','-fillpage','-append',psfile);
        end
%
% Plots Showing Differences as Color
%
        if m==1
          hf1 = figure;
          orient landscape;
          subplot(2,3,1);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,2);
          imshowpair(vp(:,:,2),reg(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,3);
          imshowpair(vp(:,:,2),tf1);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
        else
          figure(hf1);
          subplot(2,3,4);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,5);
          imshowpair(vp(:,:,2),reg(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,6);
          imshowpair(vp(:,:,2),tf2);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          sgtitle({['Series ' snt]; sltt{o}}, ...
                  'FontSize',16,'FontWeight','bold');
%
          print('-dpsc2','-r600','-fillpage','-append',psfile);
%
        end
%
     end
%
     v(:,:,:,l) = reg;  % Replace moving image with registered image
%
   end
%
% Save Registered Image Data to a MAT File
%
   st = stxtr{k};       % Series description

%
   matfile = ['T1rho_S' snt '.mat'];
%
   save(matfile,'ddir','fnams','iszs','nfile','nsls','nslt', ...
        'pspcs','splt','scmx','sns','snt','st','v');
%
end
%
% Write Registration Transformations to a Spreadsheet
%
dirstr = split(pwd,filesep);
dirstr = dirstr{end};
dirstr = strrep(dirstr,' ','_');
%
xlsnam = [dirstr '.xlsx'];
%
colnam1 = {'Series','MRI_Type','Registered_MRI_Times'};
%
t1 = table(string(sr),fitr,string(sltt),'VariableNames',colnam1);
%
colnam2 = {'tx_3D','ty_3D','tz_3D','rx_3D','ry_3D','rz_3D', ...
           'tx_2Dslice1','ty_2Dslice1','tx_2Dslice2','ty_2Dslice2'};
%
t2 = array2table(tmr,'VariableNames',colnam2);
%
tt = [t1 t2];
%
writetable(tt,xlsnam,'WriteMode','replacefile');
%
% Close All Windows
%
close all;
%
% Use Series Descriptions to Find T2* Series
%
ids = contains(stxt,'*')|contains(stxt,'T2','IgnoreCase',true); % T2*
%
ids = [ids; false];
idds = diff(ids);
idss = find(idds==1)+1;
idse = find(idds==-1);
%
ns = length(idss);
nse = length(idse);
if ns~=nse
  error([' *** ERROR in rd_m_dicom:  Number of starting and ending', ...
         ' indices for sets of T2* series are not equal!']);
end
%
ids = [idss idse];      % Index to start of sets (1st column) and end (2nd column)
idv = diff(ids,1,2)>2;  % Must be at least three echo times
ids = ids(idv,:);       % Valid sets of echo times
%
ns = size(ids,1);       % Number of sets of T2* sequences
%
netn = diff(ids,1,2)+1; % Number of echo times
%
clear idds idse idss nse;
%
% Check for Duplicate Echo Times
%
netmin = min(netn);
%
if netmin==max(netn)
  idd = 0;
else
  idd = find(netn>netmin);             % Index to duplicates
end
%
% Get Series Index and Echo Times as Matrices
%
idse = zeros(netmin,ns);               % Series index matrix
etnm = zeros(netmin,ns);               % Echo time into matrix
%
for k = 1:ns
%
   idsx = (ids(k,1):ids(k,2))';        % Index to series
   e = cell2mat(etn(idsx));            % Echo times for series
%
   if any(k==idd)                      % Correct duplicates
     [e,ide] = unique(e,'last');       % Get last duplicates
     if size(e,1)~=netmin
       error([' *** rd_m_dicom:  Number of echo times are not the', ...
              ' same!']);
     end
     idsx = idsx(ide);  % Correct series index
   end
%
   [~,idsrt] = sort(e);
   if all(idsrt~=(1:netmin)')
     error(' *** rd_m_dicom:  Echo times are not in order!');
   end
%
   if any(k==idd)                      % Warn about duplicates
     fprintf(1,'\n');
     warning(' *** WARNING:  Duplicate echo times!');
     fprintf(1,' Using series:\n');
     sstr = cellstr([int2str(sn(idsx)) repmat(' - ',netmin,1) ...
                     char(stxt{idsx})]);
     fprintf(1,'   Series %s\n',sstr{:})
     fprintf(1,'\n');
   end
%
   idse(:,k) = idsx;    % Series index into matrix
   etnm(:,k) = e;       % Echo time into matrix
%
end
%
netn = netmin;          % Number of echo times
%
% Get Index to Echo Time = 5 ms
%
d = etnm-5;
[~,id5] = min(d.*d);
id5 = id5';
%
% Get Additional T2* Variables
%
nfiles = nimages(ids(:,1));  % Numbers of T2* files
iszss = isz(ids(:,1),:);     % Image sizes in pixels
pspcss = pspc(ids(:,1),:);   % Pixel sizes
%
ddirss = ddirs(idse);        % Subdirectories for T2* series
afiless = afiles(idse);      % T2* files
snss = sn(idse);             % Series numbers
stxts = stxt(idse);          % T2* series
%
% Loop through the T2* Series
%
nreg = ns*(netn-1);     % Number of volume registrations
tms = zeros(nreg,10);   % 3D and 2D registration translations/rotations
ett = cell(nreg,1);     % Registration echo times
ss = cell(nreg,1);      % Series
fits = string(repmat('T2star',nreg,1));% Type of fit (T1rho or T2star)
%
for k = 1:ns
%
% Get Series Specific Indices
%
   ddir = ddirss(:,k);  % Subdirectories for series
%
   nfile = nfiles(k);   % Number of image files (slices) in this series
   etns = etnm(:,k);    % T2* echo times
   iszs = iszss(k,:);   % Image size
   pspcs = pspcss(k,:); % Pixel size
%
   sns = snss(:,k);     % Series numbers
   snt = int2str(sns(1));              % First series number as text
   psfile = ['S' snt '.ps'];           % PS file name
%
% Get Images and Maximum Scaled Image Values
%
   valmx = zeros(nfile*netn,1);
   v = zeros([iszs nfile netn]);
%
   for l = 1:netn
%
% Get Image File Names
%
      fnams = afiless{l,k};            % T2* files
%
      for m = 1:nfile
         info = dicominfo(fullfile(ddir{l},fnams{m}));
         sl = double(info.RescaleSlope);
         y0 = double(info.RescaleIntercept);
         img = dicomread(info);
         img = sl*double(img)+y0;
         v(:,:,m,l) = img;
         idx = l*nfile+m-nfile;
         valmx(idx) = max(img(:));     % Slice maximum
      end
   end
%
   scmx = 10*fix(max(valmx)/10);       % Round maximum value down
%
% Register the Different Echo Time Images to the 5 ms Echo Time Image
%
   rsl(2,1) = floor(nfile/4);
   rsl = [rsl(2,1); 3*rsl(2,1)];       % Slices for 2D registration
%
   pxmx = 3*iszs(1);                   % Maximum X pixels
   px = 0:20:pxmx;                     % Grid of 20 pixels
   nx = size(px,2);                    % Number of X grid lines
   pymx = iszs(2);                     % Maximum Y pixels
   py = 0:20:pymx;                     % Grid of 20 pixels
   ny = size(py,2);                    % Number of Y grid lines
%
   t2 = cell(2,1);      % 2D translations
%
% Register Echo Times Before 5 ms Echo Times
%
   idsrb = id5(k):-1:1;      % Echo times before 5 ms echo time
   nb = size(idsrb,2)-1;     % Number of before echo times to register
%
   for l = 1:nb
%
      n1 = idsrb(l);    % Index to fixed image
      n2 = idsrb(l+1);  % Index to moving image
%
      o = netn*k-netn-k+l+1;
      ss{o} = ['S' snt];
      ett{o} = ['Echo time ' sprintf('%.2f',etns(n2)), ...
                ' to ' sprintf('%.2f',etns(n1)) ' ms'];
%
      [reg,t] = elastix(v(:,:,:,n2),v(:,:,:,n1),[], ...
                               'Parameters_RigidBody.txt');     % 3D
%
      t = t.TransformParameters{1};

      rx = t.TransformParameters(1);
      ry = t.TransformParameters(2);
      rz = t.TransformParameters(3);
      tx = t.TransformParameters(4);
      ty = t.TransformParameters(5);
      tz = t.TransformParameters(6);
      tn = [tx; ty; tz; rad2deg*[rx; ry; rz]];
      t = sprintf(['tx = %.1f, ty = %.1f, tz = %.1f, rx =  %.1f, ', ...
                   'ry =  %.1f, rz =  %.1f'],tn);
%
      [tf1,dr1,dc1] = dftreg(v(:,:,rsl(1),n1), ...
                             v(:,:,rsl(1),n2),100);        % 2D
      t2{1} = sprintf('tx = %.1f, ty = %.1f',dc1,dr1);
      [tf2,dr2,dc2] = dftreg(v(:,:,rsl(2),n1), ...
                             v(:,:,rsl(2),n2),100);        % 2D
      t2{2} = sprintf('tx = %.1f, ty = %.1f',dc2,dr2);
%
      tms(o,:) = [tn',dc1,dr1,dc2,dr2];
%
% Plot Registration
%
      for m = 1:2
%
% Plots with Grids
%
        vp = squeeze(v(:,:,rsl(m),n1:-1:n2));
        irng = [min(vp(:)) max(vp(:))];
        figure;
        subplot(2,1,1);
        montage({vp(:,:,1),vp(:,:,2),reg(:,:,rsl(m))}, ...
                'Size',[1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        orient landscape;
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t,'FontSize',12,'FontWeight','bold');
        ylabel('3D','FontSize',12,'FontWeight','bold');
        title({['Series ' snt]; ett{o}; ['Slice ' int2str(rsl(m))]}, ...
              'FontSize',16,'FontWeight','bold');
%
        subplot(2,1,2);
        if m==1
          montage({vp(:,:,1),vp(:,:,2),tf1},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        else
          montage({vp(:,:,1),vp(:,:,2),tf2},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        end
        colormap gray;
        n3 = nb-l+1;
        brighten(n3*n3/32);
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t2{m},'FontSize',12,'FontWeight','bold');
        ylabel('2D','FontSize',12,'FontWeight','bold');
%
        if l==1&&m==1
          print('-dpsc2','-r600','-fillpage',psfile);
        else
          print('-dpsc2','-r600','-fillpage','-append',psfile);
        end
%
% Plots Showing Differences as Color
%
        if m==1
          hf1 = figure;
          orient landscape;
          subplot(2,3,1);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,2);
          imshowpair(vp(:,:,2),reg(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,3);
          imshowpair(vp(:,:,2),tf1);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
        else
          figure(hf1);
          subplot(2,3,4);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,5);
          imshowpair(vp(:,:,2),reg(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,6);
          imshowpair(vp(:,:,2),tf2);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          sgtitle({['Series ' snt]; ett{o}}, ...
                  'FontSize',16,'FontWeight','bold');
%
          print('-dpsc2','-r600','-fillpage','-append',psfile);
%
        end
%
     end
%
     v(:,:,:,n2) = reg;  % Replace moving image with registered image
%
   end
%
% Register Echo Times After 5 ms Echo Times
%
   idsra = id5(k):1:netn;    % Echo times after 5 ms echo time
   na = size(idsra,2)-1;     % Number of after echo times to register
%
   for l = 1:na
%
      n1 = idsra(l);    % Index to fixed image
      n2 = idsra(l+1);  % Index to moving image
%
      o = netn*k-netn-k+l+nb+1;
      ss{o} = ['S' snt];
      ett{o} = ['Echo time ' sprintf('%.2f',etns(n2)), ...
                ' to ' sprintf('%.2f',etns(n1)) ' ms'];
%
      [reg,t] = elastix(v(:,:,:,n2),v(:,:,:,n1),[], ...
                               'Parameters_RigidBody.txt');     % 3D
%
      t = t.TransformParameters{1};

      rx = t.TransformParameters(1);
      ry = t.TransformParameters(2);
      rz = t.TransformParameters(3);
      tx = t.TransformParameters(4);
      ty = t.TransformParameters(5);
      tz = t.TransformParameters(6);
      tn = [tx; ty; tz; rad2deg*[rx; ry; rz]];
      t = sprintf(['tx = %.1f, ty = %.1f, tz = %.1f, rx =  %.1f, ', ...
                   'ry =  %.1f, rz =  %.1f'],tn);
%
      [tf1,dr1,dc1] = dftreg(v(:,:,rsl(1),n1), ...
                             v(:,:,rsl(1),n2),100);        % 2D
      t2{1} = sprintf('tx = %.1f, ty = %.1f',dc1,dr1);
      [tf2,dr2,dc2] = dftreg(v(:,:,rsl(2),n1), ...
                             v(:,:,rsl(2),n2),100);        % 2D
      t2{2} = sprintf('tx = %.1f, ty = %.1f',dc2,dr2);
%
      tms(o,:) = [tn',dc1,dr1,dc2,dr2];
%
% Plot Registration
%
      for m = 1:2
%
% Plots with Grids
%
        vp = squeeze(v(:,:,rsl(m),n1:n2));
        irng = [min(vp(:)) max(vp(:))];
        figure;
        subplot(2,1,1);
        montage({vp(:,:,1),vp(:,:,2),reg(:,:,rsl(m))}, ...
                'Size',[1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        orient landscape;
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t,'FontSize',12,'FontWeight','bold');
        ylabel('3D','FontSize',12,'FontWeight','bold');
        title({['Series ' snt]; ett{o}; ['Slice ' int2str(rsl(m))]}, ...
              'FontSize',16,'FontWeight','bold');
%
        subplot(2,1,2);
        if m==1
          montage({vp(:,:,1),vp(:,:,2),tf1},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        else
          montage({vp(:,:,1),vp(:,:,2),tf2},'Size', ...
                  [1 3],'DisplayRange',irng,'ThumbnailSize',iszs);
        end
        colormap gray;
        n3 = l+2;
        brighten(n3*n3/32);
        hold on;
        plot3(repmat(px,2,1),repmat([0; pymx],1,nx),ones(2,nx),'r-');
        plot3(repmat([0; pxmx],1,ny),repmat(py,2,1),ones(2,ny),'r-');
        xlabel(t2{m},'FontSize',12,'FontWeight','bold');
        ylabel('2D','FontSize',12,'FontWeight','bold');
%
        print('-dpsc2','-r600','-fillpage','-append',psfile);
%
% Plots Showing Differences as Color
%
        if m==1
          hf1 = figure;
          orient landscape;
          subplot(2,3,1);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,2);
          imshowpair(vp(:,:,2),reg(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,3);
          imshowpair(vp(:,:,2),tf1);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
        else
          figure(hf1);
          subplot(2,3,4);
          imshowpair(vp(:,:,1),vp(:,:,2));
          title({'Unregistered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,5);
          imshowpair(vp(:,:,2),reg(:,:,rsl(m)));
          xlabel(t,'FontSize',11,'FontWeight','bold');
          title({'3D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          subplot(2,3,6);
          imshowpair(vp(:,:,2),tf2);
          xlabel(t2{m},'FontSize',11,'FontWeight','bold');
          title({'2D registered',['Slice ', int2str(rsl(m))]}, ...
                'FontSize',12,'FontWeight','bold');
%
          sgtitle({['Series ' snt]; ett{o}}, ...
                  'FontSize',16,'FontWeight','bold');
%
          print('-dpsc2','-r600','-fillpage','-append',psfile);
%
        end
%
      end
%
      v(:,:,:,n2) = reg;  % Replace moving image with registered image
%
   end
%
% Save Registered Image Data to a MAT File
%
   fnams = afiless{:,k};               % T2* files
   nsls = nfile;                       % Number of slices
   nfile = nfile*netn;                 % Total number of files
   st = stxts{1,k};                    % Series description
%
   matfile = ['T2star_S' snt '.mat'];
%
   save(matfile,'ddir','etns','fnams','id5','iszs','nfile','nsls', ...
        'netn','pspcs','scmx','sns','snt','st','v');
end
%
% Write Registration Transformations to the Spreadsheet
%
t1 = table(string(ss),fits,string(ett),'VariableNames',colnam1);
%
t2 = array2table(tms,'VariableNames',colnam2);
%
tt = [t1 t2];
%
writetable(tt,xlsnam,'WriteVariableNames',false,'WriteMode','append');
%
return