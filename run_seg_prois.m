%
% Get Subject Directories
%
rdirs = dir('0*');
rdirs = {rdirs([rdirs.isdir]').name}';
nsubj = size(rdirs,1);
%
% nn = [1:15 17]';
% rdirs = rdirs(nn);
% nsubj = size(rdirs,1);
%
tim = zeros(nsubj,1);   % Time to create masks
%
for kk = 1:nsubj
% for kk = 2:12
% for kk = 10:nsubj
% for kk = [16 18:20]
% for kk = 19:20
% for kk = [10 13:20]
% for kk = [4 6:9 11 12]
% for kk = 11:12
%
   cd(rdirs{kk});
%
   tstart = tic;
   pwd
   seg_prois;
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