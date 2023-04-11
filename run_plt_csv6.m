%
% Get Subject Directories
%
% subjn = [3;7;12;17;18;19;20;21;24;30;31;32;34;35;38;46];
subjn = [42;48;51;60];
nsubj = size(subjn,1);
%
sdirs = int2str(subjn);
sdirs = [repmat('0',nsubj,1) sdirs];
sdirs = strrep(cellstr(sdirs),' ','0');
%
tim = zeros(nsubj,1);   % Time to perform registration
%
for kk = 1:nsubj
%
   cd(sdirs{kk});
%
   tstart = tic;
   pwd
   plt_csv6;
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