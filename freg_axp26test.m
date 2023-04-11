%#######################################################################
%
%             * Femur REGions using AXial Plane Program *
%
%          M-File which reads the femur "standard" grid and divides the
%     the femur into six regions based on an axial plane view of the 
%     trochlea.  The regions are plotted and the points in each region 
%     are written to the MS-Excel spreadsheet:  Partitions6_fe.xlsx.
%
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
%     11-Nov-2014 * Mack Gardner-Morse
%

%#######################################################################
%
% Plots?
%
iplt = true;
% iplt = false;
iplt2 = true;           % Debugging plots
% iplt2 = false;          % Debugging plots
% iprt = true;
iprt = false;
%
% Convert Radians to Degrees and Coordinate Tolerance
%
rad2deg = 180/pi;
tol = 1e-10;            % Tolerance on grid coordinates
%
% MS-Excel Spreadsheet Name
%
xnam = 'Partitions6_fe_test.xlsx';
%
% Get Cylinder Radius Mean Value
%
load fscale3.mat;
%
rcm = mode(rsf.*rc);    % Use mode to avoid round off errors
%
% Read Radii Data from MS-Excel Spreadsheet f14e.xlsx
%
[num,txt,raw] = xlsread('f14e.xls');
%
rc = num(:,7);          % Cylinder radii
%
% Data Directory and File Names
%
load fscale3e.mat;
fdir = 'Biomarker Female Femur MAT files';
d = dir(fullfile(fdir,'*EAbc7t.mat'));
fnams = deblank(char(d.name));
nf = size(fnams,1);   % Number of files
%
tdir = '..\Biomarkers Tibia Bone Mat Files\';
tnams = dir([tdir '*CS.mat']);
tnams = char(tnams.name);
%
% Get Tibia Widths
%
nft = size(tnams,1);
tpw = zeros(nft,1);
for k = 1:nft
   fnamt = tnams(k,:);
   load(fullfile(tdir,fnamt));
   tpw(k) = max(xyzpt(:,2))-min(xyzpt(:,2));
end
tnams = tnams(:,1:5);
tpw_mn = mean(tpw);
%
% Get Standard Grid
%
load fb_rnga.mat;
if iplt2
  hf0 = figure;
end
% h = quadsurf(quad,[tq zq zeros(nq,1)]);
% hold on;
%
% Posterior and Medial/Lateral Cutoffs
%
tmin = -145;            % Angle (theta) cutoff (-150, -145, or -140)
y0 = 0;                 % Y cutoff (-1, 0, or 1)
%
% Get Indices to the Posterior and Medial/Lateral Regions
%
ipos = find(tq<=tmin);  % Posterior
idl = find(zq>=y0);     % Lateral
idm = find(zq<y0);      % Medial
%
% Loop through MAT Files
%
for kf = 1:nf           % Loop through MAT files
%
   fnam = fnams(kf,:);                 % File name
   idot = strfind(fnam,'.');
   fstr = fnam(1:idot-1);              % File name w/o extension
%
% Load Data
%
%    d = dir(fullfile(ddir,[fstr(1:8) '*_sc.mat']));    % Check for modified data
%    fn = sortrows(char(d.name));
%    nfsz = size(fn,1);
%    fn = fn(nfsz,:);     % Use the most recent modified data
%
   load(fullfile(fdir,fnam));
%
% Combine Bone Data and Transform Coordinates to Femur Coordinate System
%
   datb = comb_dat(datlb,datmb,dattb);
   [datb,datcp] = coord_tf(xyzo,xyzax',datb,datcp);
%
% Reverse X-Axis for Right Knees
%
   [nsl,ns,is] = sl_info(datb);
   if leg
     for ks = 1:nsl
        datb{ks}(:,1) = -datb{ks}(:,1);
     end
     datcp{1}(:,1) = -datcp{1}(:,1);
   end
%
% Get Maximum "Peaks" of Femur in X-Direction
%
   mx = zeros(nsl,1);
   idx = zeros(nsl,1);
   my = zeros(nsl,1);
   for ks = 1:nsl
      [mx(ks) idx(ks)] = max(datb{ks}(:,1));
      my(ks) = datb{ks}(idx(ks),2);
   end
%
   [mxs ids] = sort(mx);
   idp = find(my(ids)>0);
   idp = ids(idp(end-2:end));          % Medial peak index to top three
   idn = find(my(ids)<0);
   idn = ids(idn(end-2:end));          % Lateral peak index to top three
   mxp = mean(mx(idp));                % Medial peak X
   myp = mean(my(idp));                % Medial peak Y
   mxn = mean(mx(idn));                % Lateral peak X
   myn = mean(my(idn));                % Lateral peak Y
%
% Get Trochlea Groove
%
   idxmn = [round(mean(idp)); round(mean(idn))];
   idxmn = sort(idxmn);
   idmn = find(ids>idxmn(1)&ids<idxmn(2));
   idmn = ids(idmn(1:3));
   mxm = mean(mx(idmn));
   mym = mean(my(idmn));
%
% Get Line Directions
%
   xyc = [mxm mym];
   xym = [mxp myp];
   xyl = [mxn myn];
   vm = xym-xyc;
   vl = xyl-xyc;
   if iplt2
     figure(hf0);
     plt_datsl(datb,'k.-',0.5);
     view(-90,90);
     hold on;
     plot(mx(idp),my(idp),'yo','MarkerSize',7,'LineWidth',0.5);
     plot(mx(idn),my(idn),'yo','MarkerSize',7,'LineWidth',0.5);
     plot(mx(idmn),my(idmn),'yo','MarkerSize',7,'LineWidth',0.5);
   end
%
% Get Center of Trochlea
%
%    xyc = datcp{1}(1,1:2);              % Center of the trochlea
   xyc = [0 0];         % Origin
%
% Get Trochlea Points for Each Slice
%
   ip = NaN(nsl,2);
   idx = [];
   %
   for ks = 1:nsl
      xys = mean(datb{ks}(:,1:2));
      mxs = [mx(ks) my(ks)];
      vx = 4*(xys-mxs);
      yrng = sort([my(ks); xys(2)+vx(2)]);
      t = 1.2*yrng(2)./vm(2);
      if t>=0
        ip(ks,:) = lsect3(xyc(1,1:2),t*vm,mxs,vx)';
        if iplt2
          figure(hf0);
          h1 = plot([xyc(1); xyc(1)+t*vm(1)],[xyc(2); xyc(2)+t*vm(2)],'r-','LineW',2);
          hold on;
          h2 = plot([mxs(1); mxs(1)+4*vx(1)],[mxs(2); mxs(2)+4*vx(2)],'g-','LineW',2);
          pause;
          delete([h1,h2]);
        end
      end
      if isnan(ip(ks,1))
        t = 1.2*yrng(1)./vl(2);
        if t>=0
          ip(ks,:) = lsect3(xyc(1,1:2),t*vl,mxs,vx)';
          if iplt2
            figure(hf0);
            h1 = plot([xyc(1); xyc(1)+t*vl(1)],[xyc(2); xyc(2)+t*vl(2)],'r-','LineW',2);
            hold on;
            h2 = plot([mxs(1); mxs(1)+4*vx(1)],[mxs(2); mxs(2)+4*vx(2)],'g-','LineW',2);
            pause;
            delete([h1,h2]);
          end
        end
      end
      if ~isnan(ip(ks,1))
        id = find(datb{ks}(:,1)>sign(mx(ks))*ip(ks,1));
        ids = (is(ks)+1:is(ks+1))';
        id = ids(id);
        idx = [idx; id];
      end
   end
%
% Get Femur Mesh
%
   trib = mk_tri4f(datb);
   xyzb = cell2mat(datb);
%
% Plot Femur in Cartesian Coordinates
%
   if iplt
     hfc = figure;
     orient landscape;
     view(3);
     hold on;
     hm = trimesh(trib,xyzb(:,1),xyzb(:,2),xyzb(:,3),'FaceColor', ...
                  'none','EdgeColor','b','LineWidth',0.5);
     ht = plot3(xyzb(idx,1),xyzb(idx,2),xyzb(idx,3),'ro', ...
                'LineWidth',0.5);
     axis equal;
     xlabel({'X (mm)';['\leftarrow posterior / anterior', ...
            ' \rightarrow']},'FontSize',12,'FontWeight','bold');
     ylabel({'Y (mm)';'\leftarrow lateral / medial \rightarrow'}, ...
            'FontSize',12,'FontWeight','bold');
     zlabel({'Z (mm)';'\leftarrow inferior / superior \rightarrow'}, ...
            'FontSize',12,'FontWeight','bold');
     title(fstr,'Interpreter','none','FontSize',16,'FontWeight','bold');
     if iprt
       print('-dpsc2','-r300',[fstr(1:5) '.ps']);
     end
   end
%
% Convert to Polar Coordinates
%
   datbp = xzpl2pol(datb);
   trzb = cell2mat(datbp);
   trzb(:,1) = trzb(:,1)*rad2deg;
%
% Scale Theta and "Z" Coordinates
%
   tsc = rc(kf)./rcm;
%    tqs = tsc*tq;        % Scale by radii
   tqs = tq;            % No scaling of theta
   iz = strmatch(fstr(1:5),tnams);
   zsc = tpw(iz)./tpw_mn;
   zqs = zsc*zq;
%
% Find Boundary Points on Trochlea
%
   ir = nod2tri(idx,trib,2);
   ina = in_tri2d(trib,trzb(:,[1 3]),[tqs zqs]);
   ina = find(ina);
   in = in_tri2d(trib(ir,:),trzb(:,[1 3]),[tqs zqs]);
   in = find(in);
   itta = nod2tri(ina,trig,2);
   itt = nod2tri(in,trig,2);
   bida = meshbnd4(trig(itta,:));
   bid = meshbnd4(trig(itt,:));
   idb = setdiff(bid,bida);
   idbl = intersect(idl,idb);          % Grid points in lateral compartment
   idbm = intersect(idm,idb);          % Grid points in medial compartment
%
% Get Unique Z Values in the Lateral and Medial Compartments
%
   zs = sort(zqs);
   izu = find(diff(zs)>tol);
   zu = zs([izu; izu(end)+1]);         % Unique Z values
   nu = size(zu,1);
%
   izl = find(zu>=y0);
   zul = zu(izl);
   nul = size(zul,1);
%
   izm = find(zu<y0);
   zum = zu(izm);
   num = size(zum,1);
%
% Extend Lines Based on Boundary Points
%
   vl = [zqs(idbl) ones(size(idbl))];
   lpl = vl\tqs(idbl);                 % Slope and intercept of theta as a function of Z
   imxl = find(zu>max(zqs(idbl)));
   tcol = [tqs(idbl); [zu(imxl) ones(size(imxl))]*lpl];    % Theta cutoffs
   zcol = [zqs(idbl); zu(imxl)];
%
   vm = [zqs(idbm) ones(size(idbm))];
   lpm = vm\tqs(idbm);                 % Slope and intercept of theta as a function of Z
   imxm = find(zu<min(zqs(idbm)));
   tcom = [tqs(idbm); [zu(imxm) ones(size(imxm))]*lpm];    % Theta cutoffs
   zcom = [zqs(idbm); zu(imxm)];
%
% Find Grid Points with Theta > Theta Cut Offs
%
   ztl = sortrows([zcol tcol]);
   [zlu izlu] = unique(ztl(:,1));
   ztl = ztl(izlu,:);
   idtl = [];
   for kl = 1:nul
      zuc = zul(kl);
      idlc = find(ztl(:,1)>zuc-tol&ztl(:,1)<zuc+tol);
      if isempty(idlc)
        tuc = [zuc 1]*lpl;
      else
        tuc = ztl(idlc,2);
        tuc = tuc(1);
      end
      idtl = [idtl; find(zqs>zuc-tol&zqs<zuc+tol&tqs>=tuc)];
   end
%
   ztm = sortrows([zcom tcom]);
   [zmu izmu] = unique(ztm(:,1));
   ztm = ztm(izmu,:);
   idtm = [];
   for km = 1:num
      zuc = zum(km);
      idmc = find(ztm(:,1)>zuc-tol&ztm(:,1)<zuc+tol);
      if isempty(idmc)
        tuc = [zuc 1]*lpm;
      else
        tuc = ztm(idmc,2);
        tuc = tuc(1);
      end
      idtm = [idtm; find(zqs>zuc-tol&zqs<zuc+tol&tqs>=tuc)];
   end
%
% Plot Femur in Polar Coordinates
%
   if iplt
     hfp = figure;
     orient landscape;
     view(2);
     hold on;
     hm = trimesh(trib,trzb(:,1),trzb(:,3),'Color','b', ...
                  'LineWidth',0.5);
     hg = plot(tqs,zqs,'k.','MarkerSize',6,'LineWidth',0.5);
     hql = plot(tqs(idtl),zqs(idtl),'ro','MarkerSize',4, ...
                'LineWidth',0.5);
     hqm = plot(tqs(idtm),zqs(idtm),'go','MarkerSize',4, ...
                'LineWidth',0.5);
   end
%
% Get Remaining Indices to the Six (6) Regions
%
   iant = [idtl; idtm]; % Anterior (Trochlea)
   ictr = true(nq,1);
   ictr(iant) = false;
   ictr(ipos) = false;
   ictr = find(ictr);   % Center
%
   idl_pos = intersect(idl,ipos);      % Lateral posterior region
   idl_ctr = intersect(idl,ictr);      % Lateral center region
   idl_ant = idtl;                     % Lateral anterior region
%
   idm_pos = intersect(idm,ipos);      % Medial posterior region
   idm_ctr = intersect(idm,ictr);      % Medial center region
   idm_ant = idtm;                     % Medial anterior region
%
% Plot Regions
%
   if iplt
     figure(hfp);
%      plot(tqs(idl_ant),zqs(idl_ant),'ro');
%      plot(tqs(idm_ant),zqs(idm_ant),'mo');
%      plot(tqs(idl_ctr),zqs(idl_ctr),'gs');
%      plot(tqs(idm_ctr),zqs(idm_ctr),'bs');
     plot(tqs(idl_pos),zqs(idl_pos),'y^','MarkerSize',4, ...
          'LineWidth',0.5);
     plot(tqs(idm_pos),zqs(idm_pos),'c^','MarkerSize',4, ...
          'LineWidth',0.5);
     axis tight;
     xlabel({'\Theta (degrees)';['\leftarrow posterior / anterior', ...
            ' \rightarrow']},'FontSize',12,'FontWeight','bold');
     ylabel({'Z (mm)';'\leftarrow medial / lateral \rightarrow'}, ...
            'FontSize',12,'FontWeight','bold');
     title(fstr,'Interpreter','none','FontSize',16,'FontWeight','bold');
     if iprt
       print('-dpsc2','-r300','-append',[fstr(1:5) '.ps']);
     end
   end
%
% Create Index for the Six Regions (1-6)
%
   id = ones(nq,1);
   id(idl_pos) = 1;
   id(idl_ctr) = 2;
   id(idl_ant) = 3;
   id(idm_pos) = 4;
   id(idm_ctr) = 5;
   id(idm_ant) = 6;
%
% Write Index to Spreadsheet
%
   lbl = mat2cell([repmat('Pt ',nq,1) int2str((1:nq)')],ones(nq,1));
   xlswrite(xnam,lbl,fstr(1:5),'A1');
   xlswrite(xnam,id,fstr(1:5),'B1');
%
end
%
return