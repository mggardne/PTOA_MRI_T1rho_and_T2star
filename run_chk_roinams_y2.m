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
          '016'
          '017'
          '019'
          '020'
          '021'
          '023'
          '024'
          '027'
          '028'
          '031'
          '033'
          '034'
          '035'
          '036'
          '037'
          '038'
          '040'
          '042'
          '043'
          '045'
          '046'
          '047'
          '048'
          '049'
          '050'
          '051'
          '057'
          '060'
          '061'
          '063'
          '064'
          '068'
          '069'
          '070'
          '071'
          '072'
          '074'
          '075'
          '076'
          '077'
          '079'
          '080'
          '081'
          '082'
          '083'
          '084'
          '085'
          '086'
          '087'
          '088'
          '091'
          '092'
          '093'
          '094'
          '095'
          '096'
          '097'
          '099' };
nsubj = size(sdirs,1);
bdirs = {'BASELINE';'YEAR1';'YEAR2'};
ndir = size(bdirs,1);
%
tim = zeros(nsubj,1);   % Time to check ROI names
%
for kk = 1:nsubj
%
%    for ll = 1:ndir
   for ll = 3           % YEAR2 only
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