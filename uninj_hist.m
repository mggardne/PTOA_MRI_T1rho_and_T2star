%#######################################################################
%
%  * UNINJured Knees of PTOA Year 2 Subjects qMRI HISTograms Program *
%
%          M-File which reads the qMRI results data and separates the
%     injured and uninjured knees results.  Histograms of the T1rho and
%     T2* uninjured knees results are plotted and saved to PDF files:
%     T1rho_histogram.pdf and T2star_histogram.pdf.  Histograms of the
%     T1rho and T2* injured knees results are plotted and saved to PDF
%     files: T1rho_injured_histogram.pdf and
%     T2star_injured_histogram.pdf.
%
%     NOTES:  1.  qMRI results MAT file must must be in the subdirectory
%             results\mri_fitpsvy2\ in the file:  mri_fitpsvy2.mat.
%
%             2.  PTOA injured leg data is in the CSV file:
%             results\mri_fitpsvy2\PTOA_Subjects_InjuredLeg.csv.
%
%     12-May-2026 * Mack Gardner-Morse
%

%#######################################################################
%
% Get qMRI Results
%
load results\mri_fitpsvy2\mri_fitpsvy2.mat;
%
% Get Index to Uninjured Knees in PTOA Year 2 Visit
%
iy2 = ismember(tinj.subj,y2);          % Year 2 subjects
uninj_y2 = uninj(iy2);  % Uninjured year 2 knees (0 - left, 1 - right)
inj_y2 = double(~logical(uninj_y2));   % Injured year 2 knees (0 - left, 1 - right)
idu = uninj_y2+1;       % Index to uninjured knees (1 - left, 2 - right)
idi = inj_y2+1;         % Index to injured knees (1 - left, 2 - right)
%
% T1rho Results for Year 2 Visit Uninjured Knees
%
t1r_uninj = squeeze(t1r_respx(:,3,idu));
t1r_uninj = t1r_uninj(:);
n = cellfun('length',t1r_uninj);
n = sum(n)
t1r_uninj_all = cell2mat(t1r_uninj);
n = size(t1r_uninj_all,1)
max(t1r_uninj_all)
min(t1r_uninj_all)
t1r_uninj_all = t1r_uninj_all(t1r_uninj_all<=100);
nr = size(t1r_uninj_all,1)
figure;
histogram(t1r_uninj_all);
title('PTOA Year 2 Visit Uninjured Leg T1\rho','FontSize',16, ...
      'FontWeight','bold');
xlabel('T1\rho (ms)','FontSize',16,'FontWeight','bold');
ylabel('N','FontSize',16,'FontWeight','bold');
orient landscape;
ht = text(65,80000,['Number of pixels = ',int2str(nr)], ...
          'FontSize',12,'FontWeight','bold');
print('-dpdf','-r600','-fillpage','T1rho_histogram.pdf');
%
% T2* Results for Year 2 Visit Uninjured Knees
%
t2s_uninj = squeeze(t2s_respx(:,3,idu));
t2s_uninj = t2s_uninj(:);
n = cellfun('length',t2s_uninj);
n = sum(n)
t2s_uninj_all = cell2mat(t2s_uninj);
n = size(t2s_uninj_all,1)
max(t2s_uninj_all)
min(t2s_uninj_all)
t2s_uninj_all = t2s_uninj_all(t2s_uninj_all<=65);
ns = size(t2s_uninj_all,1)
figure;
histogram(t2s_uninj_all);
title('PTOA Year 2 Visit Uninjured Leg T2*','FontSize',16, ...
      'FontWeight','bold');
xlabel('T2* (ms)','FontSize',16,'FontWeight','bold');
ylabel('N','FontSize',16,'FontWeight','bold');
orient landscape;
ht = text(45,75000,['Number of pixels = ',int2str(ns)], ...
          'FontSize',12,'FontWeight','bold');
print('-dpdf','-r600','-fillpage','T2star_histogram.pdf');
%
% T1rho Results for Year 2 Visit Injured Knees
%
t1r_inj = squeeze(t1r_respx(:,3,idi));
t1r_inj = t1r_inj(:);
n = cellfun('length',t1r_inj);
n = sum(n)
t1r_inj_all = cell2mat(t1r_inj);
n = size(t1r_inj_all,1)
max(t1r_inj_all)
min(t1r_inj_all)
t1r_inj_all = t1r_inj_all(t1r_inj_all<=100);
nr = size(t1r_inj_all,1)
figure;
histogram(t1r_inj_all);
title('PTOA Year 2 Visit Injured Leg T1\rho','FontSize',16, ...
      'FontWeight','bold');
xlabel('T1\rho (ms)','FontSize',16,'FontWeight','bold');
ylabel('N','FontSize',16,'FontWeight','bold');
orient landscape;
ht = text(65,80000,['Number of pixels = ',int2str(nr)], ...
          'FontSize',12,'FontWeight','bold');
print('-dpdf','-r600','-fillpage','T1rho_injured_histogram.pdf');
%
% T2* Results for Year 2 Visit Injured Knees
%
t2s_inj = squeeze(t2s_respx(:,3,idi));
t2s_inj = t2s_inj(:);
n = cellfun('length',t2s_inj);
n = sum(n)
t2s_inj_all = cell2mat(t2s_inj);
n = size(t2s_inj_all,1)
max(t2s_inj_all)
min(t2s_inj_all)
t2s_inj_all = t2s_inj_all(t2s_inj_all<=65);
ns = size(t2s_inj_all,1)
figure;
histogram(t2s_inj_all);
title('PTOA Year 2 Visit Injured Leg T2*','FontSize',16, ...
      'FontWeight','bold');
xlabel('T2* (ms)','FontSize',16,'FontWeight','bold');
ylabel('N','FontSize',16,'FontWeight','bold');
orient landscape;
ht = text(45,75000,['Number of pixels = ',int2str(ns)], ...
          'FontSize',12,'FontWeight','bold');
print('-dpdf','-r600','-fillpage','T2star_injured_histogram.pdf');