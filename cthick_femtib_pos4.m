%#######################################################################
%
%        * Cartilage THiCKness FEMur TIBia POSition Program 4 *
%
%          M-File which loads the PTOA ROI files to compare the
%     relative position/orientation of the tibia relative to the femur
%     and femur relative to the tibia at the time of MRI data
%     acquisition.
%
%          Plots of the femurs and tibiae are in cthick_femtib_pos4_*.p*
%     files, where * is the PTOA subject identifier and the two
%     character visit label.  The resulting relative positions are in
%     cthick_tibfem_pos4_*.xlsx (tibia with respect to the femur) and
%     cthick_femtib_pos4_*.xlsx (femur with respect to the tibia) files,
%     where the * is the visit label (BL = Baseline, Y1 = Year 1, and
%     Y2 = Year 2).
%
%     NOTES:  1.  Program must start in the PTOA main directory with
%             subject subdirectories starting with "0".
%
%             2.  PTOA ROI MAT files *_prois.mat must be in the visit
%             subdirectories under the PTOA subject directories.  T1rho
%             MAT files must start with "T1rho_S" and T2* MAT files
%             must start with "T2star_S".  See seg_prois_b.m.
%
%             3.  The M-files coord_tf.m and plt_datsl.m must be in the
%             current path or directory.
%
%             4.  Based on cthick_femtib_pos.m by Daniel R. Sturnick on
%             21 January 2014 and cthick_femtib_pos3.m.
%
%       03-Nov-2025 * Mack Gardner-Morse
%

%#######################################################################
%
% Convert Radians to Degrees
%
ang90 = pi/2;           % 90 degrees in radians
rad2deg = 180/pi;       % Conversion from radians to degrees
%
% Base Name for Postscript Plot Files
%
resdir = fullfile('Results','FemTib_Pos');
psnam = 'cthick_femtib_pos4_';         % Base Postscript file name
psnam = fullfile(resdir,psnam);        % Add results directory
pstyp = '.ps';          % File extension (type)
%
% Get Subject and Visit Directories
%
vdirs = {'BASELINE';'YEAR1';'YEAR2'};  % Visit subdirectories
vlbls = ['BL'; 'Y1'; 'Y2'];
ivisit = menu('Pick visit subdirectory',vdirs);
vdir = vdirs{ivisit};   % Visit directory
vlbl = vlbls(ivisit,:); % Visit two character label
visit = ivisit-1;       % Visit identifier - 0 = BASELINE, 1 = YEAR1, and 2 = YEAR2
%
if ivisit==1
  sdirs = dir('0*');
  sdirs = {sdirs([sdirs.isdir]').name}';
elseif ivisit==2
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
            '096' };    % Year 1 PTOA subjects
else
  sdirs = { '005'
            '006'
            '007'
            '008'
            '010'       % Right leg T1ρ/T2* scan had a quality/bloom issue in medial compartment and are unusable
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
            '033'       % T1ρ right leg scan had a quality/bloom issue and is unusable
            '034'
            '035'
            '036'
            '037'       % T2* right leg has excessive motion artifact and is unusable
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
            '076'       % T1ρ right leg scan is missing
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
            '099' };    % Year 2 PTOA subjects
end
% 
nsubj = size(sdirs,1);
%
% Get Indices for Sex (0 - Female and 1 - Male) and Injured Leg
% (0 - Left and 1 - Right)
%
xlsnam = 'PTOA Cohort - Sex and Injured Limb.xlsx';
t = readtable(xlsnam);
%
ptoa_subj = t.ptoa_subj;
ptoa_subj = strrep(ptoa_subj,'PTOA','');
ptoa_subj = str2double(ptoa_subj);
%
sex = 2-t.dem01;
injleg = 2-t.dem34;
%
% Initialize Figures
%
hf1 = figure;
hf2 = figure;
hf3 = figure;
%
% Tibial Axes' Origin Markers and Axis Vector Colors
%
%   First Column = T1rho (ires = 0)
%   Second Column = T2* (ires = 1)
%   First Row = Injured Leg (ileg = injleg)
%   Second Row = Uninjured Leg (ileg ~= injleg)
%
mrks = {'gs' 'ro'; 'm^' 'k*'};
qcols = {'r' 'g'; 'k' 'm'};
%
% Total Number of ROI Files
%
nf = 4*nsubj;           % T1rho and T2* for both left and right knees
%
% Labels for Legs
%
tlegs = {'Left'; 'Right'};
%
% Initialize Output Variables
%
ap = NaN(nf,1);
ml = NaN(nf,1);
si = NaN(nf,1);
ir = NaN(nf,1);
vv = NaN(nf,1);
fe = NaN(nf,1);
%
ap2 = NaN(nf,1);
ml2 = NaN(nf,1);
si2 = NaN(nf,1);
ir2 = NaN(nf,1);
vv2 = NaN(nf,1);
fe2 = NaN(nf,1);
%
% Analysis Identifiers
%     Column 1 = PTOA subject number
%     Column 2 = Visit (0 = Baseline, 1 = Year 1, and 2 = Year 2)
%     Column 3 = Sex (0 = female and 1 = male)
%     Column 4 = Injured Leg (0 = left and 1 = right)
%     Column 5 = Analysis (0 = T1rho and 1 = T2*)
%     Column 6 = Leg (0 = left and 1 = right)
%
ids = NaN(nf,6);
%
% Loop through Subjects
%
for ks = 1:nsubj        % Loop through PTOA subjects
%
   sdir = sdirs{ks};
   subjn = str2double(sdir);           % Subject number
%
   psnamo = [psnam sdir '_' vlbl pstyp];    % Postscript output file name
%
   idx = find(subjn==ptoa_subj);       % Get subject identifiers
   injlegi = injleg(idx);              % Get subject injured leg
   id = [subjn visit sex(idx) injlegi];
%
   svdir = fullfile(sdir,vdir);        % Subject with visit directory
%
% Get ROI Segmentation MAT File Names
%
   d = dir(fullfile(svdir,'T*01_prois.mat'));
   roinams = {d.name}';                % 3D ROI MAT files
   nroi = size(roinams,1);
%
% Loop through MAT Files
%
   for km = 1:nroi      % Loop through MAT files
%
      roinam = deblank(roinams{km,:}); % File name
      ires = str2double(roinam(2))-1;  % ires = 0 - T1rho, ires = 1 - T2*
      idu = strfind(roinam,'_');
      ttxt = roinam(1:idu(1)-1);       % Text for plot titles
%
% Load ROI MAT File
%
      load(fullfile(svdir,roinam),'f3','fs','ibone','leg','t3', ...
           'rmat','rtmat','xyzc','xyzt');
%
      ttxt = [ttxt ' Series ' fs(2:end)];        % Text for plot titles
%
% Get Leg
%
      ileg = double(strcmp('R',leg));  % ileg = 0 - left, ileg = 1 - right
      ttxt = [ttxt ' ' tlegs{ileg+1} ' Leg'];    % Text for plot titles
%
% Get Index and Analysis Identifiers
%
      kf = 4*ks+2*ires+ileg-3;         % Index into results variables
      ids(kf,:) = [id ires ileg];      % Get analysis identifiers
%
% Get Center of Femoral Coordinate System as Vectors
%
      xcf = repmat(xyzc(1),3,1);
      ycf = repmat(xyzc(2),3,1);
      zcf = repmat(xyzc(3),3,1);
%
%  Get Center of Tibial Coordinate System as Vectors
%
      xct = repmat(xyzt(1),3,1);
      yct = repmat(xyzt(2),3,1);
      zct = repmat(xyzt(3),3,1);
%
% Plot Cordinate System Axes and Data in mRI Scanner Coordinate System
%
      figure(hf1);
      orient landscape;
      view(3);
      hold on;
      axis equal;
%
      datb = f3(2,:)';  % Femur bone segmentation slices
      plt_datsl(datb(ibone(:,1)),'b');
      axis equal;
      plot3(xyzc(:,1),xyzc(:,2),xyzc(:,3),'go','MarkerSize',7, ...
            'MarkerFaceColor','g');
      quiver3(xcf,ycf,zcf,rmat(1,:)',rmat(2,:)',rmat(3,:)',50, ...
              'Color',[0 0.6 0],'LineWidth',2);
      quiver3(xcf,ycf,zcf,rmat(1,:)',rmat(2,:)',rmat(3,:)',50, ...
              'Color',[0 0.6 0],'LineWidth',2);
%
      dattb = t3(2,:)'; % Tibia bone segmentation slices
      plt_datsl(dattb(ibone(:,2)),'m');
      plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'ro','MarkerSize',7, ...
            'MarkerFaceColor','r');
      quiver3(xct,yct,zct,rtmat(1,:)',rtmat(2,:)',rtmat(3,:)',50, ...
              'Color',[0.6 0 0],'LineWidth',2);
      quiver3(xct,yct,zct,rtmat(1,:)',rtmat(2,:)',rtmat(3,:)',50, ...
              'Color',[0.6 0 0],'LineWidth',2);
%
      title({ttxt; 'Position of Tibia Relative to Femur'}, ...
            'FontSize',20,'FontWeight','bold','Interpreter','none');
      xlabel('X (Scanner)','FontSize',16,'FontWeight','bold');
      ylabel('Y (Scanner)','FontSize',16,'FontWeight','bold');
      zlabel('Z (Scanner)','FontSize',16,'FontWeight','bold');
      if km==1&&ks==1
        print('-dpsc2','-r600','-fillpage',psnamo);
      else
        print('-dpsc2','-r600','-fillpage','-append',psnamo);
      end
%
% Transform Tibial Coordinate System Into Femoral Coordinate System
%
      [xyzt_trf,xyzf_trf] = coord_tf(xyzc,rmat',{xyzt},{xyzc});
      xyzt_trf = xyzt_trf{:};
      xyzf_trf = xyzf_trf{:};
%   
      xcf = repmat(xyzf_trf(1),3,1);
      ycf = repmat(xyzf_trf(2),3,1);
      zcf = repmat(xyzf_trf(3),3,1);
%  
      xct = repmat(xyzt_trf(1),3,1);
      yct = repmat(xyzt_trf(2),3,1);
      zct = repmat(xyzt_trf(3),3,1);
%
      xyzaxt_trf = zeros(3);
      xyzaxt_trf(1,:) = rtmat(1,:)*rmat';
      xyzaxt_trf(1,:) = xyzaxt_trf(1,:)./norm(xyzaxt_trf(1,:));
      xyzaxt_trf(2,:) = rtmat(2,:)*rmat';
      xyzaxt_trf(2,:) = xyzaxt_trf(2,:)./norm(xyzaxt_trf(2,:));
      xyzaxt_trf(3,:) = rtmat(3,:)*rmat';
      xyzaxt_trf(3,:) = xyzaxt_trf(3,:)./norm(xyzaxt_trf(3,:));
%
      xyzaxf_trf = [1 0 0; 0 1 0; 0 0 1];   % Transformed femur rotation matrix
%
% Plot Both Coordinate Systems Transformed To the Femoral Coordinates
%
      figure(hf2);
%
% Plot Femur
%
      if km==1
        orient landscape;
        hold on;
        plot3(xcf,ycf,zcf,'ks','MarkerSize',7,'LineWidth',1);
        quiver3(xcf,ycf,zcf,xyzaxf_trf(:,1),xyzaxf_trf(:,2), ...
                xyzaxf_trf(:,3),50,'b','LineWidth',3);
      end
%
% Plot Tibiae
%
      irow = 2-(ileg==injlegi);
      mrk = mrks{irow,ires+1};
      qcol = qcols{irow,ires+1};
      plot3(xct,yct,zct,mrk,'MarkerSize',7,'LineWidth',1);
      quiver3(xct,yct,zct,xyzaxt_trf(:,1),xyzaxt_trf(:,2), ...
              xyzaxt_trf(:,3),50,'Color',qcol,'LineWidth',2);
%
% Add Title and Axis Labels
%
      if km==nroi
        view(145,30);
        axis equal;
        title({'Position of Tibia Relative to Femur';['Injured ',  ...
              'T1rho (Green&Red)/Uninjured T1rho (Red&Green)']; ...
              ['Injured T2* (Magenta&Black)/Uninjured T2* ', ...
              '(Black&Magenta)']},'FontSize',16,'FontWeight','bold', ...
              'Interpreter','none');
        xlabel('X (Posterior-Anterior)', 'FontSize',12,'FontWeight', ...
               'bold','Interpreter','none');
        ylabel('Y (Medial-Lateral)', 'FontSize',12,'FontWeight', ...
               'bold','Interpreter','none');
        zlabel('Z (Inferior-Superior)', 'FontSize',12,'FontWeight', ...
               'bold','Interpreter','none');
      end
%
%  Save Position of the Tibia Relative to the Femur
%
      ap(kf) = xyzt_trf(1);
      ml(kf) = xyzt_trf(2);
      si(kf) = xyzt_trf(3);
%   
%  Save Orientation of the Tibia Relative to the Femur   
%  
      xaxt = xyzaxt_trf(:,1);
      yaxt = xyzaxt_trf(:,2);
      zaxt = xyzaxt_trf(:,3);
      irot = atan2(yaxt(2),yaxt(1))-ang90;
      ir(kf) = irot*rad2deg;
      flex = ang90-atan2(zaxt(3),zaxt(1));
      fe(kf) = flex*rad2deg;
      vval = atan2(zaxt(3),zaxt(2))-ang90;
      vv(kf) = vval*rad2deg;   
%
% Transform Femoral Coordinate System Into the Tibial Coordinate System
%
      [xyzf_trt,xyzt_trt] = coord_tf(xyzt,rtmat',{xyzc},{xyzt});
      xyzf_trt = xyzf_trt{:};
      xyzt_trt = xyzt_trt{:};
%   
      xcf = repmat(xyzf_trt(1),3,1);
      ycf = repmat(xyzf_trt(2),3,1);
      zcf = repmat(xyzf_trt(3),3,1);
%  
      xct = repmat(xyzt_trt(1),3,1);
      yct = repmat(xyzt_trt(2),3,1);
      zct = repmat(xyzt_trt(3),3,1);
%
      xyzaxf_trt = zeros(3);
      xyzaxf_trt(1,:) = rmat(1,:)*rtmat';
      xyzaxf_trt(1,:) = xyzaxf_trt(1,:)./norm(xyzaxf_trt(1,:));
      xyzaxf_trt(2,:) = rmat(2,:)*rtmat';
      xyzaxf_trt(2,:) = xyzaxf_trt(2,:)./norm(xyzaxf_trt(2,:));
      xyzaxf_trt(3,:) = rmat(3,:)*rtmat';
      xyzaxf_trt(3,:) = xyzaxf_trt(3,:)./norm(xyzaxf_trt(3,:));
%
      xyzaxt_trt = [1 0 0; 0 1 0; 0 0 1];   % Transformed tibia rotation matrix
%
% Plot Both Coordinate Systems Transformed to the Tibial Coordinates
%
   figure(hf3);
%
% Plot Tibia
%
      if km==1
        orient landscape;
        hold on;
        plot3(xct,yct,zct,'ro','MarkerSize',6,'LineWidth',1);
        quiver3(xct,yct,zct,xyzaxt_trt(:,1),xyzaxt_trt(:,2), ...
                xyzaxt_trt(:,3),50,'k','LineWidth',3);
      end
%
% Plot Femurs
%
      irow = 2-(ileg==injlegi);
      mrk = mrks{irow,ires+1};
      qcol = qcols{irow,ires+1};
      plot3(xcf,ycf,zcf,mrk,'MarkerSize',6,'LineWidth',1);
      quiver3(xcf,ycf,zcf,xyzaxf_trt(:,1),xyzaxf_trt(:,2), ...
              xyzaxf_trt(:,3),50,'Color',qcol,'LineWidth',2);
%
% Add Title and Axis Labels
%
      if km==nroi
        view(145,30);
        axis equal;
        title({'Position of Femur Relative to Tibia';['Injured ',  ...
              'T1rho (Green&Red)/Uninjured T1rho (Red&Green)']; ...
              ['Injured T2* (Magenta&Black)/Uninjured T2* ', ...
              '(Black&Magenta)']},'FontSize',16,'FontWeight','bold', ...
              'Interpreter','none');
        xlabel('X (Posterior-Anterior)', 'FontSize',12,'FontWeight', ...
               'bold','Interpreter','none');
        ylabel('Y (Medial-Lateral)', 'FontSize',12,'FontWeight', ...
               'bold','Interpreter','none');
        zlabel('Z (Inferior-Superior)', 'FontSize',12,'FontWeight', ...
               'bold','Interpreter','none');
      end
%
%  Save Position of the Femur Relative to the Tibia
%
      ap2(kf) = xyzf_trt(1);
      ml2(kf) = xyzf_trt(2);
      si2(kf) = xyzf_trt(3);
%   
%  Save Orientation of the Tibia Relative to the Femur   
%  
      xaxf = xyzaxf_trt(:,1);
      yaxf = xyzaxf_trt(:,2);
      zaxf = xyzaxf_trt(:,3);
      irot = atan2(yaxf(2),yaxf(1))-ang90;
      ir2(kf) = irot*rad2deg;
      flex = ang90-atan2(zaxf(3),zaxf(1));
      fe2(kf) = flex*rad2deg;
      vval = atan2(zaxf(3),zaxf(2))-ang90;
      vv2(kf) = vval*rad2deg;   
%
% Clear Figure 1
%
      clf(hf1);
%   
   end
%
% Save Plots in Figures 2 and 3
%
   figure(hf2);
   print('-dpsc2','-r600','-fillpage','-append',psnamo);
%
   figure(hf3);
   print('-dpsc2','-r600','-fillpage','-append',psnamo);
%
% Clear Figures 2 and 3
%
   clf(hf2);
   clf(hf3);
%
end
%
% Close Plot Windows
%
close all;
%
% Put Data for Relative Position of the Tibia w.r.t. Femur into Tables
%
hdrs = {'Subject','Visit','Sex','InjLeg','Result','Leg'};
t0 = array2table(ids,'VariableNames',hdrs);
t1 = table(ap,ml,si,ir,fe,vv);
t = [t0 t1];            % Table of relative positions
%
% Breakdown by Results and Injured Leg
%
idt1r = t.Result==0;    % T1rho results
tt1r = t(idt1r,:);      % T1rho results
idt2s = t.Result==1;    % T2* results
tt2s = t(idt2s,:);      % T2* results
%
% T1rho Injured - Uninjured Legs
%
idi = tt1r.InjLeg==tt1r.Leg;           % Injured legs
sids = tt1r.Subject(idi);              % Subjects with an injured leg
ninj = size(sids,1);                   % Number of injured legs
%
idl = tt1r.Subject==repmat(sids,1,size(tt1r,1))';     % Index to both legs
%
idu = xor(idl&repmat(idi,1,ninj),idl); % Uninjured legs
idu = any(idu')';                      % Uninjured legs
%
hdrs_diff = {'ap_diff','ml_diff','si_diff','ir_diff','fe_diff', ...
             'vv_diff'};
idi2 = idi;
if ivisit==3
  idi2(91) = false;     % Subject 76 right leg missing T1ρ scan
end
tdr = array2table(tt1r{idi2,7:12}-tt1r{idu,7:12},'VariableNames', ...
                  hdrs_diff);
%
tab_hdrs = tt1r.Properties.VariableNames;
lenhdrs = cellfun('length',tab_hdrs);
hdrsi = replaceBetween(tab_hdrs,lenhdrs,lenhdrs+1,'_Inj', ...
                       'Boundaries','exclusive');
hdrsu = replaceBetween(tab_hdrs,lenhdrs,lenhdrs+1,'_Uninj', ...
                       'Boundaries','exclusive');
%
idd = [1:3,5];
hdrsd = hdrs(idd);
lenhdrsd = cellfun('length',hdrsd);
hdrsd = replaceBetween(hdrsd,lenhdrsd,lenhdrsd+1,'_Diff', ...
                       'Boundaries','exclusive');
%
tdr = [renamevars(tt1r(idu,idd),hdrs(idd),hdrsd) tdr];     % Add identifiers
%
tdr = [renamevars(tt1r(idi2,:),tab_hdrs,hdrsi) ...
       renamevars(tt1r(idu,:),tab_hdrs,hdrsu) tdr];
%
% T2* Injured - Uninjured Legs
%
idi = tt2s.InjLeg==tt2s.Leg;           % Injured legs
sids = tt2s.Subject(idi);              % Subjects with an injured leg
ninj = size(sids,1);                   % Number of injured legs
%
idl = tt2s.Subject==repmat(sids,1,size(tt2s,1))';     % Index to both legs
%
idu = xor(idl&repmat(idi,1,ninj),idl); % Uninjured legs
idu = any(idu')';                      % Uninjured legs
%
tds = array2table(tt2s{idi,7:12}-tt2s{idu,7:12},'VariableNames', ...
                  hdrs_diff);
%
tds = [renamevars(tt2s(idu,idd),hdrs(idd),hdrsd) tds];     % Add identifiers
%
tds = [renamevars(tt2s(idi,:),tab_hdrs,hdrsi) ...
       renamevars(tt2s(idu,:),tab_hdrs,hdrsu) tds];
%
% Output Tables to MS-Excel Spreadsheet
%
xlsnam = ['cthick_tibfem_pos4_' vlbl '.xlsx'];
xlsnam = fullfile(resdir,xlsnam);
%
shtnams = {'Tibia wrt Femur Data'; 'T1rho Injured-Uninjured'; ...
           'T2star Injured-Uninjured'};
%
t = fillmissing(t,'constant',-1);      % -1 for missing values
%
writetable(t,xlsnam,'Sheet',shtnams{1});
writetable(tdr,xlsnam,'Sheet',shtnams{2});
writetable(tds,xlsnam,'Sheet',shtnams{3});
%
% Save Results to a MAT File
%
matnam = fullfile(resdir,['cthick_tibfem_pos4_' vlbl '.mat']);
save(matnam,'idi','idi2','ids','idu','ap','ml','si','ir','fe','vv', ...
     't','tdr','tds');
%
% Save Data for Relative Position of the Femur w.r.t. Tibia
%
t2 = table(ap2,ml2,si2,ir2,fe2,vv2);
t = [t0 t2];            % Table of relative positions
%
% Breakdown by Results and Injured Leg
%
idt1r = t.Result==0;    % T1rho results
tt1r = t(idt1r,:);      % T1rho results
idt2s = t.Result==1;    % T2* results
tt2s = t(idt2s,:);      % T2* results
%
% T1rho Injured - Uninjured Legs
%
idi = tt1r.InjLeg==tt1r.Leg;           % Injured legs
sids = tt1r.Subject(idi);              % Subjects with an injured leg
ninj = size(sids,1);                   % Number of injured legs
%
idl = tt1r.Subject==repmat(sids,1,size(tt1r,1))';     % Index to both legs
%
idu = xor(idl&repmat(idi,1,ninj),idl); % Uninjured legs
idu = any(idu')';                      % Uninjured legs
%
hdrs_diff = {'ap2_diff','ml2_diff','si2_diff','ir2_diff', ...
             'fe2_diff','vv2_diff'};
idi2 = idi;
if ivisit==3
  idi2(91) = false;     % Subject 76 right leg missing T1ρ scan
end
tdr = array2table(tt1r{idi2,7:12}-tt1r{idu,7:12},'VariableNames', ...
                  hdrs_diff);
%
tab_hdrs = tt1r.Properties.VariableNames;
lenhdrs = cellfun('length',tab_hdrs);
hdrsi = replaceBetween(tab_hdrs,lenhdrs,lenhdrs+1,'_Inj', ...
                       'Boundaries','exclusive');
hdrsu = replaceBetween(tab_hdrs,lenhdrs,lenhdrs+1,'_Uninj', ...
                       'Boundaries','exclusive');
%
tdr = [renamevars(tt1r(idu,idd),hdrs(idd),hdrsd) tdr];     % Add identifiers
%
tdr = [renamevars(tt1r(idi2,:),tab_hdrs,hdrsi) ...
       renamevars(tt1r(idu,:),tab_hdrs,hdrsu) tdr];
%
% T2* Injured - Uninjured Legs
%
idi = tt2s.InjLeg==tt2s.Leg;           % Injured legs
sids = tt2s.Subject(idi);              % Subjects with an injured leg
ninj = size(sids,1);                   % Number of injured legs
%
idl = tt2s.Subject==repmat(sids,1,size(tt2s,1))';     % Index to both legs
%
idu = xor(idl&repmat(idi,1,ninj),idl); % Uninjured legs
idu = any(idu')';                      % Uninjured legs
%
tds = array2table(tt2s{idi,7:12}-tt2s{idu,7:12},'VariableNames', ...
                  hdrs_diff);
%
tds = [renamevars(tt2s(idu,idd),hdrs(idd),hdrsd) tds];     % Add identifiers
%
tds = [renamevars(tt2s(idi,:),tab_hdrs,hdrsi) ...
       renamevars(tt2s(idu,:),tab_hdrs,hdrsu) tds];
%
% Output Tables to MS-Excel Spreadsheet
%
xlsnam = ['cthick_femtib_pos4_' vlbl '.xlsx'];
xlsnam = fullfile(resdir,xlsnam);
%
t = fillmissing(t,'constant',-1);      % -1 for missing values
%
shtnams{1} = 'Femur wrt Tibia Data';   % Femur wrt to tibia data
%
writetable(t,xlsnam,'Sheet',shtnams{1});
writetable(tdr,xlsnam,'Sheet',shtnams{2});
writetable(tds,xlsnam,'Sheet',shtnams{3});
%
% Save Results to a MAT File
%
matnam = fullfile(resdir,['cthick_femtib_pos4_' vlbl '.mat']);
save(matnam,'idi','idi2','ids','idu','ap2','ml2','si2','ir2','fe2', ...
     'vv2','t','tdr','tds');
%
return