%
% Get Subject Directories
%
sdirs = dir('0*');
sdirs = {sdirs([sdirs.isdir]').name}';
nsubj = size(sdirs,1);
%
tim = zeros(nsubj,1);   % Time to perform registration
%
% Loop through Subject Directories
%
for kk = 1:nsubj
%
   cd(sdirs{kk});
%
   tstart = tic;
   pwd
   rd_pd_dicom;         % Read proton density (PD) series
   close all;
   tim(kk) = toc(tstart);
%
   cd ..
%
end
%
tim
%
return