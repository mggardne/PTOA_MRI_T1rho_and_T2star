%
% Convert Radians to Degrees
%
rad2deg = 180/pi;
%
% Posterior and Medial/Lateral Cutoffs
%
tmin = -145/rad2deg;    % Angle (theta) cutoff (-150, -145, or -140)
y0 = 0;                 % Y cutoff (-1, 0, or 1)
hpi = pi/2;
%
% Read Segmentations
%
rdir = 'RHO';
brois = rd_prois(fullfile(rdir,'012_L_ALL_RHO_EB.csv'));
%
% Get Femur Bone Segmentations
%
frois = brois(1).rois(2).roi;          % Femur bone segmentation
datl = frois(1).data3;                 % Lateral condyle
datm = frois(2).data3;                 % Medial condyle
datt = frois(3).data3;                 % Trochlea
%
% Approximate Origin
%
xyzc = mean([ 112.8638  -70.4279    2.2751
              114.4067  -69.8591    1.8244
              114.3033  -71.9095    2.4991 ]);
%
% Transform to Femur-Like Patient Coordinate System
%
[~,ny] = plane_fit(datl{1}(:,1),datl{1}(:,2),datl{1}(:,3));
ny = -ny'./norm(ny);
nx = cross(ny,[0 0 1]);
nx = nx./norm(nx);
nz = cross(nx,ny);
nz = nz./norm(nz);
rmat = [nx; ny; nz]';
[datlf,datmf,dattf] = coord_tf(xyzc,rmat,datl,datm,datt);
%
% Get Condyle Sagittal Data
%
xyzlat = cell2mat(datlf);
xyzmed = cell2mat(datmf);
%
% Fit Cylinder to Slice Data
%
ri = 15;                % Near minimum radius
xyz1i = mean(xyzlat);
xyz2i = mean(xyzmed);
%
% Fit Condyles with a Cylinder
%
[r,xyz1,xyz2] = cyl_fit(ri,xyz1i,xyz2i,[xyzlat; xyzmed]);
xyzlny = [xyz1; xyz2];
%
% Get Coordinate System
%
yax = -diff(xyzlny);
sy = norm(yax);         % Width of medial and lateral points along Y-axis
yax = yax./sy;
zax = [0 0 1];
xax = cross(yax,zax);
xax = xax./norm(xax);
zax = cross(xax,yax);
zax = zax./norm(zax);
%
% Get Axis Origin on Cylinder Axis Line Nearest Notch Point
% (pt2line.m)
%
xyzo = pt2line(xyz1,yax,zeros(1,3));
%
xyzax = [xax; yax; zax];            % Rotation matrix
%
% Update Femur-Like Patient Coordinate System
%
[datlf,datmf,dattf] = coord_tf(xyzo,eye(3),datlf,datmf,dattf);
%
% hf1 = figure;
% hold on;
% plt_datsl(datlf,'b.-',0.5);
% plt_datsl(datmf,'g.-',0.5,[0 0.7 0]);
% plt_datsl(dattf,'k.-',0.5);
% axis equal;
% view(3);
% grid on;
% xlabel('X (mm)','FontSize',12,'FontWeight','bold');
% ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
% zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
% view(40,20);
% %
% axlim = axis;
% pxyz = [axlim(1) 0 axlim(5); axlim(2) 0 axlim(5); ...
%         axlim(2) 0 axlim(6); axlim(1) 0 axlim(6); axlim(1) 0 axlim(5)];
% hp = patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),pxyz(:,3),'FaceColor', ...
%            [0.7 0.7 0.7],'EdgeColor','r');
%
% Convert XZ to Polar Coordinates
%
% datlbp = comb_dat(datlf,datmf,dattf);  % Left femur bone segmentation in patient CS
% xyzlbp = cell2mat(datlbp);
% [thlp,rlp,zlp] = cart2pol(xyzlbp(:,1),xyzlbp(:,3),xyzlbp(:,2));
%
% Get Rotated Femur System
%
[datlf2,datmf2,dattf2] = coord_tf(zeros(1,3),xyzax',datlf,datmf,dattf);
%
hf2 = figure;
hold on;
plt_datsl(datlf2,'b.-',0.5);
plt_datsl(datmf2,'g.-',0.5,[0 0.7 0]);
plt_datsl(dattf2,'k.-',0.5);
axis equal;
view(3);
grid on;
xlabel('X (mm)','FontSize',12,'FontWeight','bold');
ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
view(40,20);
%
axlim = axis;
pxyz = [axlim(1) 0 axlim(5); axlim(2) 0 axlim(5); ...
        axlim(2) 0 axlim(6); axlim(1) 0 axlim(6); axlim(1) 0 axlim(5)];
patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),pxyz(:,3),'FaceColor', ...
      [0.7 0.7 0.7],'EdgeColor','r');
%
% Convert XZ to Polar Coordinates
%
datlbf = comb_dat(datlf2,datmf2,dattf2);
xyzlbf = cell2mat(datlbf);
[thlf,rlf,zlf] = cart2pol(xyzlbf(:,1),xyzlbf(:,3),xyzlbf(:,2));
idc = thlf>hpi;
thlf(idc) = thlf(idc)-2*pi;
figure;
plot3(thlf,zlf,rlf,'.');
axlim = axis;
view(2);
%
dr = axlim(6)-axlim(5);
rmn = mean(axlim(5:6));
axlim(5) = rmn-0.75*dr;
axlim(6) = rmn+1.25*dr;
%
ppxyz  = [tmin axlim(3) axlim(5); tmin axlim(4) axlim(5); ...
          tmin axlim(4) axlim(6); tmin axlim(3) axlim(6); ...
          tmin axlim(3) axlim(5)];
[ppxyz(:,1),ppxyz(:,3),ppxyz(:,2)] = pol2cart(ppxyz(:,1), ...
                                              ppxyz(:,3),ppxyz(:,2));
figure(hf2);
patch(ppxyz(:,1),ppxyz(:,2),ppxyz(:,3),ppxyz(:,3),'FaceColor', ...
      [0.7 0.7 0.7],'EdgeColor','r');
%
% Right Leg
%
broisr = rd_prois(fullfile(rdir,'012_R_ALL_RHO_EB.csv'));
froisr = broisr(1).rois(2).roi;
datlr = froisr(1).data3;               % Lateral condyle
datmr = froisr(2).data3;               % Medial condyle
dattr = froisr(3).data3;               % Trochlea
%
xyz0 = mean([ -109.7968  -67.9641   10.4014
              -111.3118  -68.5811   10.7347 ]);
%
% Transform to Femur-Like Coordinate System
%
[~,ny] = plane_fit(datlr{1}(:,1),datlr{1}(:,2),datlr{1}(:,3));
ny = -ny'./norm(ny);
nx = cross(ny,[0 0 1]);
nx = nx./norm(nx);
nz = cross(nx,ny);
nz = nz./norm(nz);
rmat = [nx; ny; nz]';
[datlfr,datmfr,dattfr] = coord_tf(xyz0,rmat,datlr,datmr,dattr);
%
% Get Condyle Sagittal Data
%
xyzlatr = cell2mat(datlfr);
xyzmedr = cell2mat(datmfr);
%
% Fit Cylinder to Slice Data
%
xyz1ir = mean(xyzlatr);
xyz2ir = mean(xyzmedr);
%
% Fit Condyles with a Cylinder
%
[rr,xyz1r,xyz2r] = cyl_fit(ri,xyz1ir,xyz2ir,[xyzlatr; xyzmedr]);
xyzlnyr = [xyz1r; xyz2r];
%
% Get Coordinate System
%
yaxr = -diff(xyzlnyr);
sy = norm(yaxr);        % Width of medial and lateral points along Y-axis
yaxr = yaxr./sy;
zax = [0 0 1];
xax = cross(yax,zax);
xax = xax./norm(xax);
zax = cross(xax,yax);
zax = zax./norm(zax);
%
% Get Axis Origin on Cylinder Axis Line Nearest Notch Point
% (pt2line.m)
%
xyzor = pt2line(xyz1r,yaxr,zeros(1,3));
%
xyzaxr = [xax; yaxr; zax];             % Rotation matrix
%
% Update Femur-Like Coordinate System
%
[datlfr,datmfr,dattfr] = coord_tf(xyzor,eye(3),datlfr,datmfr,dattfr);
%
% hf3 = figure;
% hold on;
% plt_datsl(datlfr,'b.-',0.5);
% plt_datsl(datmfr,'g.-',0.5,[0 0.7 0]);
% plt_datsl(dattfr,'k.-',0.5);
% axis equal;
% grid on;
% xlabel('X (mm)','FontSize',12,'FontWeight','bold');
% ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
% zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
% view(-140,20);
%
% axlim = axis;
% pxyz = [axlim(1) 0 axlim(5); axlim(2) 0 axlim(5); ...
%         axlim(2) 0 axlim(6); axlim(1) 0 axlim(6); axlim(1) 0 axlim(5)];
% hp = patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),pxyz(:,3),'FaceColor', ...
%            [0.7 0.7 0.7],'EdgeColor','r');
%
% Convert XZ to Polar Coordinates
%
% datrbp = comb_dat(datlfr,datmfr,dattfr);    % Right femur bone segmentation in patient CS
% xyzrbp = cell2mat(datrbp);
% [thrp,rrp,zrp] = cart2pol(xyzrbp(:,1),xyzrbp(:,3),xyzrbp(:,2));
%
% Get Rotated Femur System
%
[datlf2r,datmf2r,dattf2r] = coord_tf(zeros(1,3),xyzaxr',datlfr, ...
                                     datmfr,dattfr);
%
hf4 = figure;
hold on;
plt_datsl(datlf2r,'b.-',0.5);
plt_datsl(datmf2r,'g.-',0.5,[0 0.7 0]);
plt_datsl(dattf2r,'k.-',0.5);
axis equal;
grid on;
xlabel('X (mm)','FontSize',12,'FontWeight','bold');
ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
view(-140,20);
%
axlim = axis;
pxyzr = [axlim(1) 0 axlim(5); axlim(2) 0 axlim(5); ...
         axlim(2) 0 axlim(6); axlim(1) 0 axlim(6); axlim(1) 0 axlim(5)];
hp = patch(pxyzr(:,1),pxyzr(:,2),pxyzr(:,3),pxyzr(:,3),'FaceColor', ...
           [0.7 0.7 0.7],'EdgeColor','r');
%
% Convert XZ to Polar Coordinates
%
datrbf = comb_dat(datlf2r,datmf2r,dattf2r);    % Right femur bone segmentation in patient CS
xyzrbf = cell2mat(datrbf);
xyzrbf(:,1) = -xyzrbf(:,1);            % Reflect X-axis for right knees
[thrf,rrf,zrf] = cart2pol(xyzrbf(:,1),xyzrbf(:,3),xyzrbf(:,2));
idc = thrf>hpi;
thrf(idc) = thrf(idc)-2*pi;
figure;
plot3(thrf,zrf,rrf,'.');
axlim = axis;
view(2);
%
dr = axlim(6)-axlim(5);
rmn = mean(axlim(5:6));
axlim(5) = rmn-0.75*dr;
axlim(6) = rmn+1.25*dr;
%
ppxyzr  = [tmin axlim(3) axlim(5); tmin axlim(4) axlim(5); ...
           tmin axlim(4) axlim(6); tmin axlim(3) axlim(6); ...
           tmin axlim(3) axlim(5)];
[ppxyzr(:,1),ppxyzr(:,3),ppxyzr(:,2)] = pol2cart(ppxyzr(:,1), ...
                                               ppxyzr(:,3),ppxyzr(:,2));
ppxyzr(:,1) = -ppxyzr(:,1);            % Reflect X-axis for right knees
figure(hf4);
patch(ppxyzr(:,1),ppxyzr(:,2),ppxyzr(:,3),ppxyzr(:,3),'FaceColor', ...
      [0.7 0.7 0.7],'EdgeColor','r');
%
% Get Maximum "Peaks" of Femur in X-Direction
%
nsl = sl_info(datrbf);
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
idp = ids(idp(end-2:end));             % Medial peak index to top three
idn = find(my(ids)<0);
idn = ids(idn(end-2:end));             % Lateral peak index to top three
mxp = mean(mx(idp));                   % Medial peak X
myp = mean(my(idp));                   % Medial peak Y
mxn = mean(mx(idn));                   % Lateral peak X
myn = mean(my(idn));                   % Lateral peak Y
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