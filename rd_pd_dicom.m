%#######################################################################
%
%           * ReaD Proton Density DICOM Volume Data Program *
%
%          M-File which reads the meniscus DICOM images from the proton
%     density series.  The image data is put into volume matrices.  The
%     volume matrices, and image and series information are saved to
%     separate MAT files starting with "PD_S*" with "*" representing
%     the series number.
%
%     NOTES:  None.
%
%     20-Apr-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Load Image File Data from dicom_lst2.mat
%
load dicom_lst2.mat afiles ddirs isz nimages pspc sn stxt;
%
% Use Series Descriptions to Find Proton Density (PD) Series
%
idp = contains(stxt,' PD','IgnoreCase',true);    % Proton density (PD)
np = sum(idp);
%
ddirp = ddirs(idp);     % Subdirectories for PD series
nfilep = nimages(idp);  % Numbers of PD files
afilep = afiles(idp);   % PD files
iszp = isz(idp,:);      % Image sizes in pixels
snp = sn(idp);          % Series numbers
pspcp = pspc(idp,:);    % Pixel sizes
stxtp = stxt(idp);      % PD series description
%
% Loop through the Proton Density (PD) Series
%
for k = 1:np
%
% Get Series Specific Indices
%
   ddir = ddirp{k};     % Subdirectory for series
%
   nfile = nfilep(k);   % Number of image files in this series
   nsls = nfile;        % Number of slices
   iszs = iszp(k,:);    % Image size
   pspcs = pspcp(k,:);  % Pixel size
%
   sns = snp(k);        % Series number
   snt = int2str(sns);  % Series number as text
%
% Get Image File Names
%
   fnams = afilep{k};   % PD files
%
% Get Images and Maximum Scaled Image Values
%
   valmx = zeros(nfile,1);
   v = zeros([iszs nsls]);
%
   for l = 1:nsls
      info = dicominfo(fullfile(ddir,fnams{l}));
      sl = double(info.RescaleSlope);
      y0 = double(info.RescaleIntercept);
      img = dicomread(info);
      img = sl*double(img)+y0;
      v(:,:,l) = img;
      valmx(l) = max(img(:));          % Slice maximum
   end
%
   scmx = 10*fix(max(valmx)/10);       % Round maximum value down
%
% Save Image Data to a MAT File
%
   st = stxtp{k};       % Series description
%
   matfile = ['PD_S' snt '.mat'];
%
   save(matfile,'ddir','fnams','iszs','nfile','nsls','pspcs', ...
        'scmx','sns','snt','st','v');
%
end
%
return