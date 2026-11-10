%
% Left and Right Legs
%
legs = ['L'; 'R'];
%
% Directory with MAT Files as a String
%
dirstr = split(pwd,filesep);
%
% Get Subject
%
subj = dirstr{end-1};
%
% T1rho
%
rdir = 'RHO';           % Directory for T1rho segmentations
%
for k = 1:2
%
   leg = legs(k);
%
   d = dir(fullfile(rdir,[subj '_' leg '_ALL_*.csv']));
   if length(d)>1
     warning(['Too many ROI files for: ' subj ' T1rho ' leg ' Leg!']);
   end
%
end
%
% T2*
%
rdir = 'T2S';           % Directory for T1rho segmentations
%
for k = 1:2
%
   leg = legs(k);
%
   d = dir(fullfile(rdir,[subj '_' leg '_ALL_*.csv']));
   if length(d)>1
     warning(['Too many ROI files for: ' subj ' T2* ' leg ' Leg!']);
   end
%
end
