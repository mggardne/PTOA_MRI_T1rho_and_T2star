%#######################################################################
%
%        * Cartilage THiCKness FEMur TIBia POSition Program 3 *
%
%          M-File which loads the new biomarkers femur to compare the
%     relative position/orientation of the tibia relative to the femur
%     at the time of MRI data acquisition.
%
%     NOTES:  1.  Scaled MAT files *_sc.mat must be in the directory
%             "Biomarker MAT files".  The MAT files fscale3e.mat and
%             fb_rnga.mat must be in the current path or directory.
%
%             2.  The M-files comb_dat.m, comp_msh.m in_tri2d.m,
%             mk_tri4.m, mk_tri4f.m, plane_fit.m, pol2xzpl.m, sl_dir.m,
%             sl_info.m, tri_fix2.m and xzpl2pol.m must be in the
%             current path or directory.
%
%     21Jan2014 * Daniel R. Sturnick
%                 -Edited to include EAbc7cs.mat files which are data files
%                 which only contain coordiante system data.
%
%#######################################################################
%
% Plots?
%
iplt = true;
% iplt = false;
iplt2 = false;          % Debugging plots
%
% Convert Radians to Degrees
%
rad2deg = 180/pi;  
%
% Data Directory and File Names
%
load fscale3e.mat;
d1 = dir('*EAbc7t.mat');
d2 = dir('*EAbc7cs.mat');
d = [d1; d2];
fnams = deblank(char(d.name));
nf = size(fnams,1);     % Number of files
%
% Get Index to All of the Knees
%
ida = (1:nf)';          % Index to all 26 knees
ids = (1:2:nf)';        % Index to pairs of knees
nsets = size(ids,1);    % Number of pairs (sets) of knees
%
% Load File Names and Flags
%
[fnam2,flags] = bmsubj_fil(7);
fnam1 = char(fnam2(:,1));
fnam2 = char(fnam2(:,2));
%
% Get Matching File IDs
%
[idn,idnum] = get_subj(fnams(ids,:));
[idn1,idnum1] = get_subj(fnam1);
[idn2,idnum2] = get_subj(fnam2);
idxs = zeros(nsets,1);
for k = 1:nsets
   id = find(idnum1==idnum(k));
   if isempty(id)
     id = find(idnum2==idnum(k));
   end
   idxs(k) = id;
end
%
% Get Indices for Injured (0) and Uninjured (1) Knees
%
[fs,idfs] = sort([idnum1(idxs) idnum2(idxs)]');
inj = idfs(:)-1;        % Injured = 0, Uninjured = 1
icase = [find(inj) find(~inj)];
iinj = ~inj;
%
%
%
tdir = 'Biomarkers Tibia Bone Mat Files/';
%
% Get Standard Grid
%
load fb_rnga.mat;
%
% Load figures
%
hf1 = figure;
hf2 = figure;
hf3 = figure;
%
% Initialize output variables
%
fnums = zeros(nf,1);
ap = zeros(nf,1);
ml = zeros(nf,1);
si = zeros(nf,1);
ir = zeros(nf,1);
vv = zeros(nf,1);
fe = zeros(nf,1);
ap2 = zeros(nf,1);
ml2 = zeros(nf,1);
si2 = zeros(nf,1);
ir2 = zeros(nf,1);
vv2 = zeros(nf,1);
fe2 = zeros(nf,1);
%
% Loop through MAT Files
%
for kf = 1:nf           % Loop through MAT files
% for kf = [1]          % Loop through MAT files
%
   fnam = deblank(fnams(kf,:));                 % File name
   idot = strfind(fnam,'.');
   idu = strfind(fnam,'_');
   kstr = fnam(1:idu(1)+1);
   fnum = str2double(fnam(1:idu(1)-1));
   fstr = fnam(1:idot-1);              % File name w/o extension
%
% Load Femur Data
%
%   d = dir([fstr(1:8) '*t.mat']);    % Check for modified data
%   fn = sortrows(char(d.name));
%   nfsz = size(fn,1);
%   fn = fn(nfsz,:);     % Use the most recent modified data
%
   
   load(fnam);
%
% Rename Femur Cordinate System Variables
%
   xyzaxf = xyzax;
   xyzof = xyzo;
%
%  Load Tibia Data 
%
   fnamt = strrep(fnam,'EAbc7t','TibiaCS');
   fnamt = strrep(fnamt,'EAbc7cs','TibiaCS');
   fnamt = [tdir fnamt];
%
   load(fnamt);
%
%  Rename Tibia Coordinate System Variables
%
   xyzaxt = [xyzr(:,1)'; xyzr(:,2)'; xyzr(:,3)'];
   xyzot = xyzc;
   if leg
       xyzaxt(1,:) = -1*xyzaxt(1,:);
       xyzaxf(1,:) = -1*xyzaxf(1,:);  
   end
%
% Plot Cordinate System Axes and Data in Scanner Coordinate System
%
   figure(hf1);
   hold on;
   axis equal;
   xcf = repmat(xyzof(1),3,1);
   ycf = repmat(xyzof(2),3,1);
   zcf = repmat(xyzof(3),3,1);
%  
   xct = repmat(xyzot(1),3,1);
   yct = repmat(xyzot(2),3,1);
   zct = repmat(xyzot(3),3,1);
%
   plot3(xcf,ycf,zcf,'go','MarkerSize',9,'LineWidth',2);
   hf = quiver3(xcf,ycf,zcf,xyzaxf(:,1),xyzaxf(:,2),xyzaxf(:,3),50,'b');
   set(hf,'LineWidth',3);
%
   plot3(xct,yct,zct,'ro','MarkerSize',9,'LineWidth',2);
   ht = quiver3(xct,yct,zct,xyzaxt(:,1),xyzaxt(:,2),xyzaxt(:,3),50,'k');
   set(ht,'LineWidth',3);
%
   xyzpt = xyzpt*xyzaxt + repmat(xyzot,size(xyzpt,1),1);%*xyzaxt';
   plot3(xyzpt(:,1),xyzpt(:,2),xyzpt(:,3),'k.-');
   plt_datsl(datlmtb,'b-');
   view(3);
   title({'Position of Tibia Relative to Femur';kstr}, ...
       'FontSize',16,'FontWeight','bold','Interpreter','none');
   xlabel('X (Scanner)', 'FontSize',16,'FontWeight','bold','Interpreter','none');
   ylabel('Y (Scanner)', 'FontSize',16,'FontWeight','bold','Interpreter','none');
   zlabel('Z (Scanner)', 'FontSize',16,'FontWeight','bold','Interpreter','none');
   if kf==1
       print('-dpsc2','cthick_femtib_pos');
   else
       print('-dpsc2','cthick_femtib_pos','-append');
   end
%
% Transfer Tibial Coordinate System Into Femoral Coordinate System
%
   [dataxt_trf,datpt_trf,datot_trf,dataxf_trf,datof_trf] = coord_tf(xyzof,xyzaxf',{xyzaxt},{xyzpt},{xyzot},{xyzaxf},{xyzof});
%   xyzaxt_trf = dataxt_trf{:};
   xyzpt_trf = datpt_trf{:};
   xyzot_trf = datot_trf{:};
%   xyzaxf_trf = dataxf_trf{:};
   xyzof_trf = datof_trf{:};
%
   xyzaxt_trf = zeros(3);
   xyzaxt_trf(1,:) = (xyzaxt(1,:))*xyzaxf';
   xyzaxt_trf(1,:) = xyzaxt_trf(1,:)./norm(xyzaxt_trf(1,:));
   xyzaxt_trf(2,:) = (xyzaxt(2,:))*xyzaxf';
   xyzaxt_trf(2,:) = xyzaxt_trf(2,:)./norm(xyzaxt_trf(2,:));
   xyzaxt_trf(3,:) = (xyzaxt(3,:))*xyzaxf';
   xyzaxt_trf(3,:) = xyzaxt_trf(3,:)./norm(xyzaxt_trf(3,:));
%
   xyzaxf_trf = [1 0 0; 0 1 0; 0 0 1];     
%   
   xcf = repmat(xyzof_trf(1),3,1);
   ycf = repmat(xyzof_trf(2),3,1);
   zcf = repmat(xyzof_trf(3),3,1);
%  
   xct = repmat(xyzot_trf(1),3,1);
   yct = repmat(xyzot_trf(2),3,1);
   zct = repmat(xyzot_trf(3),3,1);
%
% Plot Both Coordinate Systems Transfered To Femoral Coordinates
%
   figure(hf2)
   hold on;
   plot3(xcf,ycf,zcf,'go','MarkerSize',6,'LineWidth',1);
   hf = quiver3(xcf,ycf,zcf,xyzaxf_trf(:,1),xyzaxf_trf(:,2),xyzaxf_trf(:,3),50,'b');
   set(hf,'LineWidth',3);
%

   if any(find(iinj)==kf)
    plot3(xct,yct,zct,'go','MarkerSize',6,'LineWidth',1);
    hfi = quiver3(xct,yct,zct,xyzaxt_trf(:,1),xyzaxt_trf(:,2),xyzaxt_trf(:,3),50,'k');
    set(hfi,'LineWidth',1.5);
   else
    plot3(xct,yct,zct,'yo','MarkerSize',6,'LineWidth',1);
    hfu = quiver3(xct,yct,zct,xyzaxt_trf(:,1),xyzaxt_trf(:,2),xyzaxt_trf(:,3),50,'r');
    set(hfu,'LineWidth',1.5);
   end
   axis equal;
   title({'Position of Tibia Relative to Femur';'Injured Tibia (Black&Green)/Uninjured Tibia (Red&Yellow)'}, ...
       'FontSize',16,'FontWeight','bold','Interpreter','none');
   xlabel('X (Anterior-Posterior)', 'FontSize',12,'FontWeight','bold','Interpreter','none');
   ylabel('Y (Medial-Lateral)', 'FontSize',12,'FontWeight','bold','Interpreter','none');
   zlabel('Z (Inferior-Superior)', 'FontSize',12,'FontWeight','bold','Interpreter','none');
%
%  Save Position of Tibia Relative to Femur
%
   ap(kf) = xyzot_trf(1);
   %if leg
   %    ap(kf) = -1*ap(kf);
   %end
   ml(kf) = xyzot_trf(2);
   si(kf) = xyzot_trf(3);
%   
%  Save Orientation of the Tibia Relative to Femur   
%  
   fnums(kf) = fnum;
   xaxt = xyzaxt_trf(:,1);
   %if leg
   %    xaxt = -1*xaxt;
   %end
   yaxt = xyzaxt_trf(:,2);
   zaxt = xyzaxt_trf(:,3);
   irot = atan2(yaxt(2),yaxt(1)) - atan2(1,0);
   ir(kf) = irot*rad2deg;
   flex = atan2(1,0) - atan2(zaxt(3),zaxt(1));
   fe(kf) = flex*rad2deg;
   vval = atan2(zaxt(3),zaxt(2)) - atan2(1,0);
   vv(kf) = vval*rad2deg;   
%
% Transfer Femoral Coordinate System Into Tibial Coordinate System
%
   [dataxf_trt,datof_trt,dataxt_trt,datot_trt] = coord_tf(xyzot,xyzaxt',{xyzaxf},{xyzof},{xyzaxt},{xyzot});
%   xyzaxt_trf = dataxt_trf{:};
   xyzof_trt = datof_trt{:};
%   xyzaxf_trf = dataxf_trf{:};
   xyzot_trt = datot_trt{:};
%
    xyzaxf_trt = zeros(3);
    xyzaxf_trt(1,:) = (xyzaxf(1,:))*xyzaxt';
    xyzaxf_trt(1,:) = xyzaxf_trt(1,:)./norm(xyzaxf_trt(1,:));
    xyzaxf_trt(2,:) = (xyzaxf(2,:))*xyzaxt';
    xyzaxf_trt(2,:) = xyzaxf_trt(2,:)./norm(xyzaxf_trt(2,:));
    xyzaxf_trt(3,:) = (xyzaxf(3,:))*xyzaxt';
    xyzaxf_trt(3,:) = xyzaxf_trt(3,:)./norm(xyzaxf_trt(3,:));
%
    xyzaxt_trt = [1 0 0; 0 1 0; 0 0 1];     
%   
   xcf = repmat(xyzof_trt(1),3,1);
   ycf = repmat(xyzof_trt(2),3,1);
   zcf = repmat(xyzof_trt(3),3,1);
%  
   xct = repmat(xyzot_trt(1),3,1);
   yct = repmat(xyzot_trt(2),3,1);
   zct = repmat(xyzot_trt(3),3,1);
%
% Plot Both Coordinate Systems Transfered to Tibial Coordiantes
   figure(hf3);
   hold on;
   if any(find(iinj)==kf)
    plot3(xcf,ycf,zcf,'go','MarkerSize',6,'LineWidth',1);
    hfi = quiver3(xcf,ycf,zcf,xyzaxf_trt(:,1),xyzaxf_trt(:,2),xyzaxf_trt(:,3),50,'b');
    set(hfi,'LineWidth',1.5);
   else
    plot3(xcf,ycf,zcf,'yo','MarkerSize',6,'LineWidth',1);
    hfu = quiver3(xcf,ycf,zcf,xyzaxf_trt(:,1),xyzaxf_trt(:,2),xyzaxf_trt(:,3),50,'r');
    set(hfu,'LineWidth',1.5);
   end
%
   plot3(xct,yct,zct,'ro','MarkerSize',6,'LineWidth',1);
   ht = quiver3(xct,yct,zct,xyzaxt_trt(:,1),xyzaxt_trt(:,2),xyzaxt_trt(:,3),50,'k');
   set(ht,'LineWidth',3);
   axis equal;
   xlabel('X (Anterior-Posterior)', 'FontSize',12,'FontWeight','bold','Interpreter','none');
   ylabel('Y (Medial-Lateral)', 'FontSize',12,'FontWeight','bold','Interpreter','none');
   zlabel('Z (Inferior-Superior)', 'FontSize',12,'FontWeight','bold','Interpreter','none');
   title({'Position of Femur Relative to Tibia';'Injured Femurs(Blue&Green)/Uninjured Femurs(Red&Yellow)'}, ...
       'FontSize',16,'FontWeight','bold','Interpreter','none');
%
%  Save Position of Tibia Relative to Femur
%
   ap2(kf) = xyzof_trt(1);
   %if leg
   %    ap(kf) = -1*ap(kf);
   %end
   ml2(kf) = xyzof_trt(2);
   si2(kf) = xyzof_trt(3);
%   
%  Save Orientation of the Tibia Relative to Femur   
%  
   fnums(kf) = fnum;
   xaxf = xyzaxf_trt(:,1);
   %if leg
   %    xaxt = -1*xaxt;
   %end
   yaxf = xyzaxf_trt(:,2);
   zaxf = xyzaxf_trt(:,3);
   irot = atan2(yaxf(2),yaxf(1)) - atan2(1,0);
   ir2(kf) = irot*rad2deg;
   flex = atan2(1,0) - atan2(zaxf(3),zaxf(1));
   fe2(kf) = flex*rad2deg;
   vval = atan2(zaxf(3),zaxf(2)) - atan2(1,0);
   vv2(kf) = vval*rad2deg;   
%
% Clear Figures
%
   clf(hf1);
%   
end
%
% Save Plots
%
figure(hf2);
view(145,30)
print('-dpsc2','cthick_femtib_pos2');
figure(hf3);
view(145,30)
print('-dpsc2','cthick_femtib_pos2','-append');
%
% Save Data for Relative Position of the Tibial w.r.t. Femur
%
res = {ap ml si ir fe vv};
match = 1:nsets;
csvnam1 = 'cthick_femtib_pos_data.csv';
csvnam2 = 'cthick_femtib_pos_diffs.csv';
csvnam3 = 'cthick_femtib_pos_stats.csv';
hdr_txt1 = {'fnam' 'ap' 'ml' 'si' 'ir' 'fe' 'vv'};
hdr_txt2 = {'match' 'knum_inj' 'ap_inj' 'ml_inj' 'si_inj' 'ir_inj' 'fe_inj' 'vv_inj' ...
    'knum_uninj' 'ap_uninj' 'ml_uninj' 'si_uninj' 'ir_uninj' 'fe_uninj' 'vv_uninj' ...
    'ap_diff' 'ml_diff' 'si_diff' 'ir_diff' 'fe_diff' 'vv_diff'};
hdr_txt3 = {'' 'ap_inj' 'ml_inj' 'si_inj' 'ir_inj' 'fe_inj' 'vv_inj' ...
    'ap_uninj' 'ml_uninj' 'si_uninj' 'ir_uninj' 'fe_uninj' 'vv_uninj' ...
    'ap_diff' 'ml_diff' 'si_diff' 'ir_diff' 'fe_diff' 'vv_diff'};
fid = fopen(csvnam1,'w');
fid2 = fopen(csvnam2,'w');
fid3 = fopen(csvnam3,'w');
nc = size(hdr_txt1,2);
for k = 1:nc-1
     fprintf(fid,'"%s",',hdr_txt1{k});
end
fprintf(fid,'"%s"\n',hdr_txt1{nc});
nc = size(hdr_txt2,2);
for k = 1:nc-1
     fprintf(fid2,'"%s",',hdr_txt2{k});
end
fprintf(fid2,'"%s"\n',hdr_txt2{nc});
nc = size(hdr_txt3,2);
for k = 1:nc-1
     fprintf(fid3,'"%s",',hdr_txt3{k});
end
fprintf(fid3,'"%s"\n',hdr_txt3{nc});
res1 = [fnums res{:}];
dlmwrite(csvnam1,res1,'-append');
%
resm = [res{:}];
res_f_case_inj = resm(find(iinj),:);
res_f_case_uninj = resm(find(~iinj),:);
pdiffs = {res_f_case_inj-res_f_case_uninj};
res2 = [match' fnums(find(iinj)) res_f_case_inj fnums(find(~iinj)) res_f_case_uninj pdiffs{:}];
dlmwrite(csvnam2,res2,'-append');
%
res3 = [res_f_case_inj res_f_case_uninj pdiffs{:}];
stdev = std(res3);
mx = max(res3);
mn = min(res3);
ave = mean(res3);
stats = [ave; mx; mn; stdev];
%
[h,p] = ttest(res_f_case_inj,res_f_case_uninj,0.05,'both');
%
dlmwrite(csvnam3,stats,'-append');
dlmwrite(csvnam3,p,'-append');
%
% Save Data for Relative Position of the Femur w.r.t. Tibia
%
res = {ap2 ml2 si2 ir2 fe2 vv2};
match = 1:nsets;
csvnam1 = 'cthick_femtib_pos_data2.csv';
csvnam2 = 'cthick_femtib_pos_diffs2.csv';
csvnam3 = 'cthick_femtib_pos_stats2.csv';
hdr_txt1 = {'fnam' 'ap2' 'ml2' 'si2' 'ir2' 'fe2' 'vv2'};
hdr_txt2 = {'match' 'knum_inj' 'ap_inj' 'ml_inj' 'si_inj' 'ir_inj' 'fe_inj' 'vv_inj' ...
    'knum_uninj' 'ap_uninj' 'ml_uninj' 'si_uninj' 'ir_uninj' 'fe_uninj' 'vv_uninj' ...
    'ap_diff' 'ml_diff' 'si_diff' 'ir_diff' 'fe_diff' 'vv_diff'};
hdr_txt3 = {'' 'ap_inj' 'ml_inj' 'si_inj' 'ir_inj' 'fe_inj' 'vv_inj' ...
    'ap_uninj' 'ml_uninj' 'si_uninj' 'ir_uninj' 'fe_uninj' 'vv_uninj' ...
    'ap_diff' 'ml_diff' 'si_diff' 'ir_diff' 'fe_diff' 'vv_diff'};
fid = fopen(csvnam1,'w');
fid2 = fopen(csvnam2,'w');
fid3 = fopen(csvnam3,'w');
nc = size(hdr_txt1,2);
for k = 1:nc-1
     fprintf(fid,'"%s",',hdr_txt1{k});
end
fprintf(fid,'"%s"\n',hdr_txt1{nc});
nc = size(hdr_txt2,2);
for k = 1:nc-1
     fprintf(fid2,'"%s",',hdr_txt2{k});
end
fprintf(fid2,'"%s"\n',hdr_txt2{nc});
nc = size(hdr_txt3,2);
for k = 1:nc-1
     fprintf(fid3,'"%s",',hdr_txt3{k});
end
fprintf(fid3,'"%s"\n',hdr_txt3{nc});
res1 = [fnums res{:}];
dlmwrite(csvnam1,res1,'-append');
%
resm = [res{:}];
res_f_case_inj = resm(find(iinj),:);
res_f_case_uninj = resm(find(~iinj),:);
pdiffs = {res_f_case_inj-res_f_case_uninj};
res2 = [match' fnums(find(iinj)) res_f_case_inj fnums(find(~iinj)) res_f_case_uninj pdiffs{:}];
dlmwrite(csvnam2,res2,'-append');
%
res3 = [res_f_case_inj res_f_case_uninj pdiffs{:}];
stdev = std(res3);
mx = max(res3);
mn = min(res3);
ave = mean(res3);
stats = [ave; mx; mn; stdev];
%
[h,p] = ttest(res_f_case_inj,res_f_case_uninj,0.05,'both');
%
dlmwrite(csvnam3,stats,'-append');
dlmwrite(csvnam3,p,'-append');

fclose all;
%
save cthick_femtib_pos ap ml si ir fe vv fnums;
