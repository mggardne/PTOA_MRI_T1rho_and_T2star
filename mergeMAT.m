%#######################################################################
%
%                      * MERGE MAT File Program *
%
%          M-File which reads two mri_fitps MAT files and merges the
%     data of the first five subjects in the first MAT file with the
%     last 15 subjects in the second MAT file.  The resulting data is
%     saved into the first MAT file.
%
%     NOTES:  1.  MAT files mri_fitps.mat and mri_fitps2.mat files must
%             be in the results directory:  Results\mri_fitps\
%
%     07-Apr-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Load Data
%
s1 = load('Results\mri_fitps\mri_fitps.mat');    % First five subjects
s2 = load('Results\mri_fitps\mri_fitps2.mat');   % Last 15 subjects
%     t1r_nps: {6-D cell}
%     t1r_npx: [6-D double]
%     t1r_res: [6-D double]
%   t1r_respx: {6-D cell}
%     t1r_rss: [6-D double]
%   t1r_rsspx: {6-D cell}
%     t2s_nps: {6-D cell}
%     t2s_npx: [6-D double]
%     t2s_res: [6-D double]
%   t2s_respx: {6-D cell}
%     t2s_rss: [6-D double]
%   t2s_rsspx: {6-D cell}
%
% Merge T1rho Data
%
t1r_nps = s1.t1r_nps;
t1r_nps2 = s2.t1r_nps;
t1r_nps(6:20,:,:,:,:,:) = t1r_nps2(6:20,:,:,:,:,:);
%
t1r_npx = s1.t1r_npx;
t1r_npx2 = s2.t1r_npx;
t1r_npx(6:20,:,:,:,:,:) = t1r_npx2(6:20,:,:,:,:,:);
%
t1r_res = s1.t1r_res;
t1r_res2 = s2.t1r_res;
t1r_res(6:20,:,:,:,:,:) = t1r_res2(6:20,:,:,:,:,:);
%
t1r_respx = s1.t1r_respx;
t1r_respx2 = s2.t1r_respx;
t1r_respx(6:20,:,:,:,:,:) = t1r_respx2(6:20,:,:,:,:,:);
%
t1r_rss = s1.t1r_rss;
t1r_rss2 = s2.t1r_rss;
t1r_rss(6:20,:,:,:,:,:) = t1r_rss2(6:20,:,:,:,:,:);
%
t1r_rsspx = s1.t1r_rsspx;
t1r_rsspx2 = s2.t1r_rsspx;
t1r_rsspx(6:20,:,:,:,:,:) = t1r_rsspx2(6:20,:,:,:,:,:);
%
t2s_nps = s1.t2s_nps;
t2s_nps2 = s2.t2s_nps;
t2s_nps(6:20,:,:,:,:,:) = t2s_nps2(6:20,:,:,:,:,:);
%
t2s_npx = s1.t2s_npx;
t2s_npx2 = s2.t2s_npx;
t2s_npx(6:20,:,:,:,:,:) = t2s_npx2(6:20,:,:,:,:,:);
%
t2s_res = s1.t2s_res;
t2s_res2 = s2.t2s_res;
t2s_res(6:20,:,:,:,:,:) = t2s_res2(6:20,:,:,:,:,:);
%
t2s_respx = s1.t2s_respx;
t2s_respx2 = s2.t2s_respx;
t2s_respx(6:20,:,:,:,:,:) = t2s_respx2(6:20,:,:,:,:,:);
%
t2s_rss = s1.t2s_rss;
t2s_rss2 = s2.t2s_rss;
t2s_rss(6:20,:,:,:,:,:) = t2s_rss2(6:20,:,:,:,:,:);
%
t2s_rsspx = s1.t2s_rsspx;
t2s_rsspx2 = s2.t2s_rsspx;
t2s_rsspx(6:20,:,:,:,:,:) = t2s_rsspx2(6:20,:,:,:,:,:);
%
resdir = 'Results\mri_fitps\';
%
save(fullfile(resdir,'mri_fitps.mat'),'t1r_res','t1r_npx','t1r_rss', ...
     't1r_respx','t1r_rsspx','t1r_nps','t2s_res','t2s_npx', ...
     't2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
%
return