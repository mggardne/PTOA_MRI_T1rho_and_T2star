%
% Get Subject Directories
%
rdirs = { '005'
          '006'
          '007'
          '008'
          '010'         % Right leg T1ρ/T2* scan had a quality/bloom issue in medial compartment and are unusable
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
          '033'         % T1ρ right leg scan had a quality/bloom issue and is unusable
          '034'
          '035'
          '036'
          '037'         % T2* right leg has excessive motion artifact and is unusable
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
          '076'         % T1ρ right leg scan is missing
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
nsubj = size(rdirs,1);
bdirs = {'BASELINE';'YEAR1';'YEAR2'};
ndir = size(bdirs,1);
%
% Index to Updated T2* Subjects
%
irdirs = [ 7
          18
          31
          44
          45
          46
          48
          51
          54
          66 ]';
%
irdirs = irdirs(end-3:end);  % Only last four subjects
% irdirs = irdirs(1:end-4);    % Skip subjects in common with T1rho subjects
%
tim = zeros(nsubj,1);   % Time to create masks
%
% for kk = 1:nsubj
for kk = 24:nsubj
% for kk = [55:57 61 62]
% for kk = [37 39 63 65 66]
% for kk = irdirs         % Updated T2*
% for kk = 63:nsubj       % Updated T1rho and T2*
%
%    for ll = 1:ndir
   for ll = 3           % YEAR2 only
%
      cd(fullfile(rdirs{kk},bdirs{ll}));
%
      pwd
%
      tstart = tic;
      seg_prois_b;
      close all;
      clear brois d* f* i* l* m* nm* nr* nsl* p* rdir rm* rs* rt* rx s*;
      clear t t2 t3* tp* v x*;
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