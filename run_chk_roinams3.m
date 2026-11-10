%
% Get Subject Directories
%
% subjn = [(5:8)';10;11;14;15;17;20;21;23;27];
% nsubj = size(subjn,1);
%
% sdirs = int2str(subjn);
% sdirs = [repmat('0',nsubj,1) sdirs];
% sdirs = strrep(cellstr(sdirs),' ','0');
sdirs = { '005'
          '006'
          '007'
          '008'
          '010'
          '011'
          '014'
          '015'
          '017'
          '018'
          '020'
          '021'
          '023'
          '024'
          '025'
          '027'
          '028'
          '031'
          '033'
          '034'
          '035'
          '037'
          '038'
          '040'
          '042'
          '049'
          '051'
          '053'
          '055'
          '056'
          '057'
          '058'
          '060'
          '061'
          '062'
          '063'
          '064'
          '067'
          '070'
          '077'
          '080'
          '082'
          '091'
          '092'
          '096' };      % Year 1 PTOA subjects
%
% sdirs = dir('0*')
% sdirs = {sdirs([sdirs.isdir]').name}';
%
nsubj = size(sdirs,1);
bdirs = {'BASELINE';'YEAR1';'YEAR2'};
ndir = size(bdirs,1);
%
tim = zeros(nsubj,1);   % Time to check ROI names
%
for kk = 1:nsubj
%
%    for ll = 1:ndir
   for ll = 2           % YEAR1 only
%
      cd(fullfile(sdirs{kk},bdirs{ll}));
%
      tstart = tic;
      pwd
      chk_roinams;
      fprintf(1,'\n');
      tim(kk) = toc(tstart);
%
      cd ..\..
%
   end
%
end
%
tim
%
return