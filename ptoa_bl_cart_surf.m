%
% Plot Surfaces?
%
iplt = true;            % Plots!
% iplt = false;           % No plots
%
% Leg Labels
%
llbls = {'Left'; 'Right'};             % Leg labels
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','HSS');    % Results directory
%
% Get Subject Directories
%
sdirs = { '005'
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
nsubj = size(sdirs,1);
%
% Visits
%
vdirs = {'BASELINE';'YEAR1';'YEAR2'};
vds = {'BL';'Y1';'Y2'};
nvisit = size(vdirs,1);
%
% T1rho ROI File Names
%
rdir = 'T1rho_S*01_prois.mat';
%
% Loop Through Baseline Subjects
%
for kk = 1:nsubj
%
   sdir = sdirs{kk};
   subj = str2double(sdir);
%
% Baseline Visits Only
%
%    for ll = 1:nvisit
   for ll = 1           % BASELINE only
%
      vdir = vdirs{ll};
      vd = vds{ll};
      pdir = fullfile(sdir,vdir);      % Get path to subject's baseline visit
      fprintf(1,'Working on %s (%i of 66) ...\n',pdir,kk);
%
% Get ROI MAT Files in Subject's Baseline Visit
%
      roinams = dir(fullfile(pdir,rdir));
      roinams = {roinams.name}';
      nrois = size(roinams,1);
%
% Loop Through ROI Files
%
      for mm = 1:nrois
%
% Read ROI File
%
         roinam = roinams{mm};
         load(fullfile(pdir,roinam),'brois','f3f','leg','plan_nvs', ...
              'plan_pts','rmat','rslf','rslt','rtmat','t3tt', ...
              'tplan_nvs','tplan_pts','xyzc','xyzt');
%
% Get Femoral Mesh and Coordinates
%
         trifc = mk_tri4f(f3f(1,:)');  % Cartilage connectivity
         xyzfc = cell2mat(f3f(1,:)');  % Cartilage coordinates
         trifc = tri_fix2(trifc,xyzfc);
%
         trifb = mk_tri4f(f3f(2,:)');  % Bone connectivity
         xyzfb = cell2mat(f3f(2,:)');  % Bone coordinates
         trifb = tri_fix2(trifb,xyzfb);
%
% Get Tibia Compartments
%
         idcl = brois(2).rois(1).roi(1).imageno; % Cartilage lateral compartment
         idbl = brois(2).rois(2).roi(1).imageno; % Bone lateral compartment
         idlat = intersect(idcl,idbl)';     % Tibia lateral compartment slices
         idlat = ismember(rslt,idlat);
%
% Get Tibia Mesh and Coordinates
%
         tritcl = mk_tri6(t3tt(1,idlat)');  % Lateral cartilage connectivity
         xyztcl = cell2mat(t3tt(1,idlat)'); % Lateral cartilage coordinates
         tritcl = tri_fix2(tritcl,xyztcl);
%
         tritcm = mk_tri6(t3tt(1,~idlat)');      % Medila cartilage connectivity
         xyztcm = cell2mat(t3tt(1,~idlat)');     % Medial cartilage coordinates
         tritcm = tri_fix2(tritcm,xyztcm);
%
         tritbl = mk_tri6(t3tt(2,idlat)');  % Lateral bone connectivity
         xyztbl = cell2mat(t3tt(2,idlat)'); % Lateral bone coordinates
         tritbl = tri_fix2(tritbl,xyztbl);
%
         tritbm = mk_tri6(t3tt(2,~idlat)');      % Medial bone connectivity
         xyztbm = cell2mat(t3tt(2,~idlat)');     % Medial bone coordinates
         tritbm = tri_fix2(tritbm,xyztbm);
%
% Plot Femur
%
         if iplt
           figure;
           tri_shr(trifc,xyzfc);
           hold on;
           tri_shr(trifb,xyzfb);
           view(3);
           orient landscape;
           axis equal;
           xlabel('X (mm)','FontSize',12,'FontWeight','bold');
           ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
           zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
           ttxt = {['Subject ' sdir ' Baseline T1\rho ']; ...
                   [llbls{startsWith(llbls,leg)} ' Knee']};
           title(ttxt,'FontSize',16,'FontWeight','bold');
         end
%
         if mm==1
           ressdir = fullfile(resdir,sdir);
           if ~exist(ressdir,'dir')
             mkdir(resdir,sdir);
           end
         end
%
         resnam = fullfile(ressdir,[sdir '_' vd '_T1R_' leg]);
%
         if iplt
           print([resnam '_Femur.pdf'],'-dpdf','-bestfit');
           savefig([resnam '_Femur']);
           close(gcf);
         end
%
% Plot Tibia
%
         if iplt
           figure;
           tri_shr(tritcl,xyztcl);
           hold on;
           tri_shr(tritcm,xyztcm);
           tri_shr(tritbl,xyztbl);
           tri_shr(tritbm,xyztbm);
           view(3);
           orient landscape;
           axis equal;
           xlabel('X (mm)','FontSize',12,'FontWeight','bold');
           ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
           zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
           ttxt = {['Subject ' sdir ' Baseline T1\rho ']; ...
                   [llbls{startsWith(llbls,leg)} ' Knee']};
           title(ttxt,'FontSize',16,'FontWeight','bold');
%
           print([resnam '_Tibia.pdf'],'-dpdf','-bestfit');
           savefig([resnam '_Tibia']);
           close(gcf);
         end
%
% Femur Cartilage Segmentations
%
% 1 - Lateral
% 2 - Medial
% 3 - Trochlea
%
         idlat = brois(1).rois(1).roi(1).imageno;
         idmed = brois(1).rois(1).roi(2).imageno;
         if size(brois(1).rois(1).roi,2)>2
           idtroch = brois(1).rois(1).roi(3).imageno;
         else
           idtroch = 0;
         end
         idl = ismember(rslf,idlat);
         idm = ismember(rslf,idmed);
         idt = ismember(rslf,idtroch);
%
% Femur Cartilage Segmentations Points
%
         f3fc = f3f(1,:)';
         npts = cellfun('length',f3fc);
         idregc = zeros(sum(npts),1);
%
         if idl(1)
           id1 = sum(npts(idl));
           idregc(1:id1) = 1;
           id2 = sum(npts(idt));
           idregc(id1+1:id1+id2) = 3;
           id3 = sum(npts(idm));
           idregc(id1+id2+1:id1+id2+id3) = 2;
         else
           id1 = sum(npts(idm));
           idregc(1:id1) = 2;
           id2 = sum(npts(idt));
           idregc(id1+1:id1+id2) = 3;
           id3 = sum(npts(idl));
           idregc(id1+id2+1:id1+id2+id3) = 1;
         end
%
% Femur Bone Segmentations
%
% 1 - Lateral
% 2 - Medial
% 3 - Trochlea
%
         idlat = brois(1).rois(2).roi(1).imageno;
         idmed = brois(1).rois(2).roi(2).imageno;
         if size(brois(1).rois(2).roi,2)>2
           idtroch = brois(1).rois(2).roi(3).imageno;
         else
           idtroch = 0;
         end
         idl = ismember(rslf,idlat);
         idm = ismember(rslf,idmed);
         idt = ismember(rslf,idtroch);
%
% Femur Bone Segmentations Points
%
         f3fb = f3f(2,:)';
         npts = cellfun('length',f3fb);
         idregb = zeros(sum(npts),1);
%
         if idl(1)
           id1 = sum(npts(idl));
           idregb(1:id1) = 1;
           id2 = sum(npts(idt));
           idregb(id1+1:id1+id2) = 3;
           id3 = sum(npts(idm));
           idregb(id1+id2+1:id1+id2+id3) = 2;
         else
           id1 = sum(npts(idm));
           idregb(1:id1) = 2;
           id2 = sum(npts(idt));
           idregb(id1+1:id1+id2) = 3;
           id3 = sum(npts(idl));
           idregb(id1+id2+1:id1+id2+id3) = 1;
         end
%
% Save MAT File
%
         save([resnam '.mat'],'idregb','idregc','leg','plan_nvs', ...
              'plan_pts','rmat','rtmat','sdir','subj','tplan_nvs', ...
              'tplan_pts','trifb','tritbl','tritbm','trifc', ...
              'tritcl','tritcm','vd','vdir','xyzfb','xyztbl', ...
              'xyztbm','xyzc','xyzfc','xyzt','xyztcl','xyztcm');
%
      end
%
   end
%
end