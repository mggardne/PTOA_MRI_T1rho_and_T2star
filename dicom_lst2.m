%#######################################################################
%
%                       * DICOM LiST 2 Program *
%
%          M-File which reads a series of DICOM directories and
%      collects information about the MRI series and DICOM files.
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
%             than selects the subdirectories for processing (typically
%             the "Select all" button).
%
%             5.  A table of information is displayed to the screen,
%             written to a MS-Excel spreadsheet and to a Matlab MAT
%             file.  The spreadsheet and MAT file are written in the
%             same directory as the parent directory.  The MAT file
%             contains additional variables and all the file names in
%             each series and patient position in MRI coordinates.
%
%     06-May-2021 * Mack Gardner-Morse
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
  warning(' *** WARNING in dicom_lst2:  No series selected!');
  return
end
%
ddirs = ddirs(isel);
ns = size(ddirs,1);     % Number of series
%
% Loop through the Series and Get Header Information
%
nimages = zeros(ns,1);  % Number of DICOM files
afiles = cell(ns,1);    % All DICOM file in series
afile1 = cell(ns,1);    % First DICOM files in series
afile2 = cell(ns,1);    % Last DICOM files in series
%
adur = zeros(ns,2);     % Acquisition duration (s)
adurs = cell(ns,1);     % Acquisition duration string (mm:ss)
etn = cell(ns,1);       % Echo times as numbers for UTE T2 star sequences
ets = cell(ns,1);       % Echo times as a string for UTE T2 star sequences
psz = zeros(ns,2);      % Rows and columns
im_type = cell(ns,1);   % Image type in series
isz = zeros(ns,2);      % Width and Height
sthk = zeros(ns,1);     % Slice Thickness
sspc = zeros(ns,1);     % Spacing Between Slices
pos = zeros(ns,3);      % Patient position
pspc = zeros(ns,2);     % Pixel Spacing
ptxt = cell(ns,1);      % Protocol Name
stxt = cell(ns,1);      % Series Description
sl = zeros(ns,1);       % Rescale Slope
rinterc = zeros(ns,1);  % Rescale Intercept
splt = zeros(ns,1);     % Spin lock time (Trigger Time)
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
   afile1(k) = afiles{k}(1);
   afile2(k) = afiles{k}(nimages(k));
%
   fnam = afile1{k};
   fnam = fullfile(ddir,fnam);
%
% Get DICOM Header Information from Files
%
   if exist(fnam,'file')
%
     info = dicominfo(fnam);
%
     if isfield(info,'AcquisitionDuration')
       adur(k) = info.AcquisitionDuration;
       adurs{k} = char(duration(seconds(adur(k)),'Format','mm:ss'));
     else
       adurs{k} = '00:00';
     end
     if isfield(info,'Rows')
       psz(k,:) = [info.Rows info.Columns];
       isz(k,:) = [info.Width info.Height];
     end
     if isfield(info,'SliceThickness')
       sthk(k) = info.SliceThickness;
       pspc(k,:) = info.PixelSpacing';
       sspc(k) = info.SpacingBetweenSlices;
     end
     if isfield(info,'ImageType')
       im_type{k} = info.ImageType;
     end
     if isfield(info,'ImagePositionPatient')
       pos(k,:) = info.ImagePositionPatient';
     end
     if isfield(info,'RescaleSlope')
       sl(k) = info.RescaleSlope;
       rinterc(k) = info.RescaleIntercept;
     end
     ptxt{k} = info.ProtocolName;
     stxt{k} = info.SeriesDescription;
     sn(k) = info.SeriesNumber;
     if isfield(info,'TriggerTime')
       splt(k) = info.TriggerTime;
     end
     if isfield(info,'EchoTrainLength')
       n = info.EchoTrainLength;
       if n>nimages(k)
         n = nimages(k);
       end
       et = NaN(n,1); % Echo times for this sequence
       et(1) = info.EchoTime;
       if isfield(info,'Private_2001_1025');
         ets{k} = info.Private_2001_1025;
       else
         ets{k} = num2str(et(1));
       end
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
       etn{k} = etu(idnnan);
       if size(etn{k},1)==1&&n>1
         idslash = strfind(ets{k},'/');
         if ~isempty(idslash)
           et1 = eval(ets{k}(1:idslash-1));
           det1 = eval(ets{k}(idslash+1:end));
           etn{k} = (et1:det1:(n-1)*det1+et1)';
         else
           etn{k} = eval(ets{k});
         end
       end
       ets{k} = [ets{k} '(' int2str(n) ')'];
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
idvr = contains(stxt,'TSL');
%
splcktc = cellstr(repmat('None',ns,1));     % Cell array for all series
%
if any(idvr)
  splckt = extractBetween(stxt(idvr),'TSL','ms');    % Spin lock times as text
  splckt = strrep(splckt,'_',' ');
  splckt = strtrim(splckt);
  splckt = strrep(splckt,' ',',');
  splcktc(idvr) = splckt;
end
%
% Include UTE T2star Echo Times
%
idvu = contains(splcktc,'None')&~cellfun(@isempty,ets);
%
if any(idvu)
  splcktc(idvu) = ets(idvu);
end
%
% Put Series Data into a Table and Get Column Names
%
colnams = {'Directories','Series#','SeriesDateTime', ...
           'SeriesDescription','ProtocolName','ImageType', ...
           'AcquisitionDuration','SpinLockTimes','#ofFiles', ...
           'FirstFile','LastFile','SliceThickness', ...
           'SpacingBetweenSlices','ImageRows','ImageColumns', ...
           'Width','Height','ImagePixelX','ImagePixelY', ...
           'RescaleSlope','RescaleIntercept'};
t0 = table(string(ddirs),sn,sdat,string(stxt),string(ptxt), ...
           string(im_type),string(adurs),string(splcktc),nimages, ...
           string(afile1),string(afile2),sthk,sspc,psz(:,1), ...
           psz(:,2),isz(:,1),isz(:,2),round(pspc(:,1),3), ...
           round(pspc(:,2),3),sl,rinterc,'VariableNames',colnams)
%
% Write Table to Spreadsheet
%
xlsnam = fullfile('dicom_lst2.xlsx');
writetable(t0,xlsnam);
%
% Save MAT File
%
clear ans d0 dat det1 et1 etu idnnan idslash info info1 k l m n tim;
matnam = fullfile('dicom_lst2.mat');
save(matnam);
%
return