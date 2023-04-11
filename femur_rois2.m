iplt2 = true;
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
% Read Femur Segmentations and Axis Data
%
rdir = 'RHO';
if true
brois = rd_prois(fullfile(rdir,'012_L_ALL_RHO_EB.csv'));
baxis = rd_roi6(fullfile(rdir,'012_L_AX_FEM_AB.csv'));
%
% Get Femur Bone Segmentations
%
frois = brois(1).rois(2).roi;          % Femur bone segmentation
datl = frois(1).data3;                 % Lateral condyle
datm = frois(2).data3;                 % Medial condyle
datt = frois(3).data3;                 % Trochlea
%
% Get Center Point and Proximal Femur Points
%
datcp = baxis(1).data;  % Center point
datpf = baxis(2).data;  % Proximal femur
%
% Get Femur Coordinate System and Fitted Cylinder
%
[xyzc,rmat] = f_cs_14(datl,datm,datpf,datcp,iplt2);
title({'Femur CS'; 'Left Leg'},'FontSize',16,'FontWeight','bold');
%
% Transform to Femur Coordinate System
%
[datlf,datmf,dattf] = coord_tf(xyzc,rmat,datl,datm,datt);
%
hf2 = figure;
hold on;
plt_datsl(datlf,'b.-',0.5,[0 0 0.7]);
plt_datsl(datmf,'g.-',0.5,[0 0.5 0]);
plt_datsl(dattf,'k.-',0.5);
axis equal;
view(3);
grid on;
xlabel('X (mm)','FontSize',12,'FontWeight','bold');
ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
title({'Dividing Planes in Femur CS'; 'Left Leg'},'FontSize',16, ...
      'FontWeight','bold');
view(40,20);
axis equal;
%
axlim = axis;
pxyz = [axlim(1) 0 axlim(5); axlim(2) 0 axlim(5); ...
        axlim(2) 0 axlim(6); axlim(1) 0 axlim(6); axlim(1) 0 axlim(5)];
patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),pxyz(:,3),'FaceColor', ...
      [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
%
% Convert XZ to Polar Coordinates
%
datlbf = comb_dat(datlf,datmf,dattf);
xyzlbf = cell2mat(datlbf);
%
[thlf,rlf,zlf] = cart2pol(xyzlbf(:,1),xyzlbf(:,3),xyzlbf(:,2));
idc = thlf>hpi;
thlf(idc) = thlf(idc)-2*pi;
%
figure;
plot3(thlf,zlf,rlf,'.');
axlim = axis;
view(2);
hold on;
plot([tmin tmin],axlim(3:4),'r-','LineWidth',1);
xlabel('Theta (radians)','FontSize',12,'FontWeight','bold');
ylabel('Z (mm)','FontSize',12,'FontWeight','bold');
title({'Theta-Z CS'; 'Left Leg'},'FontSize',16,'FontWeight','bold');
%
dr = axlim(6)-axlim(5);
rmn = mean(axlim(5:6));
axlim(5) = rmn-0.75*dr;
axlim(6) = rmn+1.25*dr;
%
ppxyz = [tmin axlim(3) axlim(5); tmin axlim(4) axlim(5); ...
         tmin axlim(4) axlim(6); tmin axlim(3) axlim(6); ...
         tmin axlim(3) axlim(5)];
[ppxyz(:,1),ppxyz(:,3),ppxyz(:,2)] = pol2cart(ppxyz(:,1), ...
                                              ppxyz(:,3),ppxyz(:,2));
figure(hf2);
patch(ppxyz(:,1),ppxyz(:,2),ppxyz(:,3),ppxyz(:,3),'FaceColor', ...
      [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
%
% Get Maximum "Peaks" of Femur in X-Direction
%
nsl = size(datlbf,1);
%
mx = zeros(nsl,1);
my = zeros(nsl,1);
for ks = 1:nsl
   [mx(ks),idx] = max(datlbf{ks}(:,1));
   my(ks) = datlbf{ks}(idx,2);
end
%
[~,ids] = sort(mx);
idp = find(my(ids)>0);
idp = ids(idp(end-2:end));             % Lateral peak index to top three
idn = find(my(ids)<0);
idn = ids(idn(end-2:end));             % Medial peak index to top three
mxp = mean(mx(idp));                   % Lateral peak X
myp = mean(my(idp));                   % Lateral peak Y
mxn = mean(mx(idn));                   % Medial peak X
myn = mean(my(idn));                   % Medial peak Y
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
xyc = [mxm mym];        % Center of trochlea groove at maximum X
xyl = [mxp myp];
xym = [mxn myn];
vm = xym-xyc;
vl = xyl-xyc;
%
% Get Center of Trochlea (Center of Coordinate System)
%
% xyc = [0 0 0];          % Use origin
%
% Get Lateral and Medial Planes
%
vl = [vl 0];            % Make 3D
vl = vl./norm(vl);
vm = [vm 0];            % Make 3D
vm = vm./norm(vm);
xv = [1 0 0];
%
zvl = cross(xv,vl);
vnl = cross(vl,zvl);    % Lateral plane normal
vnl = vnl./norm(vnl);
%
zvm  = cross(xv,vm);
vnm = cross(vm,zvm);    % Medial plane normal
vnm = vnm./norm(vnm);
%
% Plot Lateral and Medial Planes
%
figure(hf2);
axis equal;
axlim = axis;
%
xlmx = -vnl(2)*axlim(4)/vnl(1);
xmmx = -vnm(2)*axlim(3)/vnm(1);
%
pxyztll  = [0 0 axlim(6); 0 0 axlim(5); xlmx axlim(4) axlim(5); xlmx ...
          axlim(4) axlim(6); 0 0 axlim(6)];
%
pxyztlm  = [0 0 axlim(6); xmmx axlim(3) axlim(6); xmmx ...
          axlim(3) axlim(5); 0 0 axlim(5); 0 0 axlim(6)];
%
patch(pxyztll(:,1),pxyztll(:,2),pxyztll(:,3),pxyztll(:,3), ...
      'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
patch(pxyztlm(:,1),pxyztlm(:,2),pxyztlm(:,3),pxyztlm(:,3), ...
      'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
%
end                     % Skip left leg
%
%
%

%
% Right Leg
%
broisr = rd_prois(fullfile(rdir,'012_R_ALL_RHO_EB.csv'));
baxisr = rd_roi6(fullfile(rdir,'012_R_AX_FEM_AB.csv'));
%
% Get Femur Bone Segmentations
%
froisr = broisr(1).rois(2).roi;
datlr = froisr(1).data3;               % Lateral condyle
datmr = froisr(2).data3;               % Medial condyle
dattr = froisr(3).data3;               % Trochlea
%
% Get Center Point and Proximal Femur Points
%
datcpr = baxisr(1).data;               % Center point
datpfr = baxisr(2).data;               % Proximal femur
%
% Get Femur Coordinate System and Fitted Cylinder
%
[xyzcr,rmatr] = f_cs_14(datlr,datmr,datpfr,datcpr,iplt2);
title({'Femur CS'; 'Right Leg'},'FontSize',16,'FontWeight','bold');
%
% Transform to Femur Coordinate System
%
[datlfr,datmfr,dattfr] = coord_tf(xyzcr,rmatr,datlr,datmr,dattr);
%
% Reverse X-Axis for Right Knees
%
for ks = 1:size(datlfr,1)
   datlfr{ks}(:,1) = -datlfr{ks}(:,1);
end
%
for ks = 1:size(datmfr,1)
   datmfr{ks}(:,1) = -datmfr{ks}(:,1);
end
%
for ks = 1:size(dattfr,1)
   dattfr{ks}(:,1) = -dattfr{ks}(:,1);
end
%
hf4 = figure;
hold on;
plt_datsl(datlfr,'b.-',0.5,[0 0 0.7]);
plt_datsl(datmfr,'g.-',0.5,[0 0.5 0]);
plt_datsl(dattfr,'k.-',0.5);
axis equal;
grid on;
xlabel('X (mm)','FontSize',12,'FontWeight','bold');
ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
title({'Dividing Planes in Femur CS';'Right Leg'},'FontSize',16, ...
      'FontWeight','bold');
view(40,20);
axis equal;
%
axlim = axis;
pxyzr = [axlim(1) 0 axlim(5); axlim(2) 0 axlim(5); ...
         axlim(2) 0 axlim(6); axlim(1) 0 axlim(6); axlim(1) 0 axlim(5)];
hp = patch(pxyzr(:,1),pxyzr(:,2),pxyzr(:,3),pxyzr(:,3),'FaceColor', ...
           [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
%
% Get Plane Point and Normal
%
pp1 = [0 0 0];
pn1 = [0 1 0];
%
% Convert XZ to Polar Coordinates
%
datrbf = comb_dat(datlfr,datmfr,dattfr);    % Right femur bone segmentation in patient CS
xyzrbf = cell2mat(datrbf);
%
[thrf,rrf,zrf] = cart2pol(xyzrbf(:,1),xyzrbf(:,3),xyzrbf(:,2));
idc = thrf>hpi;
thrf(idc) = thrf(idc)-2*pi;
%
figure;
plot3(thrf,zrf,rrf,'.');
axlim = axis;
view(2);
hold on;
plot([tmin tmin],axlim(3:4),'r-','LineWidth',1);
xlabel('Theta (radians)','FontSize',12,'FontWeight','bold');
ylabel('Z (mm)','FontSize',12,'FontWeight','bold');
title({'Theta-Z CS'; 'Left Leg'},'FontSize',16,'FontWeight','bold');
%
dr = axlim(6)-axlim(5);
rmn = mean(axlim(5:6));
axlim(5) = rmn-0.75*dr;
axlim(6) = rmn+1.25*dr;
%
ppxyzr = [tmin axlim(3) axlim(5); tmin axlim(4) axlim(5); ...
          tmin axlim(4) axlim(6); tmin axlim(3) axlim(6); ...
          tmin axlim(3) axlim(5)];
[ppxyzr(:,1),ppxyzr(:,3),ppxyzr(:,2)] = pol2cart(ppxyzr(:,1), ...
                                               ppxyzr(:,3),ppxyzr(:,2));
figure(hf4);
patch(ppxyzr(:,1),ppxyzr(:,2),ppxyzr(:,3),ppxyzr(:,3),'FaceColor', ...
      [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
%
% Get Plane Point and Normal
%
pp2 = mean(ppxyzr(1:4,:));
pn2 = diff(ppxyzr(1:3,:));
pn2 = cross(-pn2(1,:),pn2(2,:));
%
% Get Maximum "Peaks" of Femur in X-Direction
%
nsl = size(datrbf,1);
%
mx = zeros(nsl,1);
my = zeros(nsl,1);
for ks = 1:nsl
   [mx(ks),idx] = max(datrbf{ks}(:,1));
   my(ks) = datrbf{ks}(idx,2);
end
%
[~,ids] = sort(mx);
idp = find(my(ids)>0);
idp = ids(idp(end-2:end));             % Lateral peak index to top three
idn = find(my(ids)<0);
idn = ids(idn(end-2:end));             % Medial peak index to top three
mxp = mean(mx(idp));                   % Lateral peak X
myp = mean(my(idp));                   % Lateral peak Y
mxn = mean(mx(idn));                   % Medial peak X
myn = mean(my(idn));                   % Medial peak Y
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
xyc = [mxm mym];        % Center of trochlea groove at maximum X
xyl = [mxp myp];
xym = [mxn myn];
vm = xym-xyc;
vl = xyl-xyc;
%
% Get Center of Trochlea (Center of Coordinate System)
%
% xyc = [0 0 0];          % Use origin
%
% Get Lateral and Medial Planes
%
vl = [vl 0];            % Make 3D
vl = vl./norm(vl);
vm = [vm 0];            % Make 3D
vm = vm./norm(vm);
xv = [1 0 0];
%
zvl = cross(xv,vl);
pnl = cross(vl,zvl);    % Lateral plane normal
pnl = pnl./norm(pnl);
ppl = [0 0 0];          % Point in lateral plane
%
zvm  = cross(xv,vm);
pnm = cross(vm,zvm);    % Medial plane normal
pnm = pnm./norm(pnm);
ppm = [0 0 0];          % Point in medial plane
%
% Plot Lateral and Medial Planes
%
figure(hf4);
axis equal;
axlim = axis;
%
xlmx = -vnl(2)*axlim(4)/vnl(1);
xmmx = -vnm(2)*axlim(3)/vnm(1);
%
pxyztl  = [0 0 axlim(6); 0 0 axlim(5); xlmx axlim(4) axlim(5); xlmx ...
          axlim(4) axlim(6); 0 0 axlim(6)];
%
pxyztm  = [0 0 axlim(6); xmmx axlim(3) axlim(6); xmmx ...
          axlim(3) axlim(5); 0 0 axlim(5); 0 0 axlim(6)];
%
patch(pxyztl(:,1),pxyztl(:,2),pxyztl(:,3),pxyztl(:,3),'FaceColor', ...
      [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
patch(pxyztm(:,1),pxyztm(:,2),pxyztm(:,3),pxyztm(:,3),'FaceColor', ...
      [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
%
return