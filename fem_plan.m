function [lmp,lmv,ppc,ppn,ppl,pnl,ppm,pnm] = fem_plan(f3,iplt,ltxt,fsav)
%FEM_PLAN  Uses femur bone sagittal plane segmentations to calculate
%          planes dividing the femur into regions of interest (ROIs).
%
%          [LMP,LMV,PPC,PPN,PPL,PNL,PPM,PNM] = FEM_PLAN(F3) Given a
%          sagittal knee segmentation cell array, F3, calculates the
%          lateral-medial dividing plane defined by a point, LMP,
%          and a normal vector, LMV, the posterior dividing plane
%          defined by a point, PPC, and a normal vector, PPN, the
%          lateral plane dividing the lateral condyle from the trochlea
%          defined by a point, PPL, and a normal vector, PNL, and a
%          medial plane dividing the medial condyle from the trochlea
%          defined by a point, PPM, and a normal vector, PNM.
%
%          [...] = FEM_PLAN(F3,IPLT,LTXT,FSAV) given a logical true,
%          IPLT, creates plots of the femoral dividing planes. LTXT is
%          a string with the name of the leg to be used in the plot
%          title.  If the file name, FSAV, is provided, the plots are
%          saved to the file, FSAV.
%
%          NOTES:  1.  The femur bone sagittal plane segmentations
%                  should be in the femoral coordinate system.
%
%                  2.  The Matlab files plnorm.m and plt_datsl.m must
%                  be in the current directory or path.
%
%                  3.  The origin, [0 0 0], is a point in the lateral-
%                  medial plane and the lateral and medial trochlear
%                  planes.
%
%          16-Dec-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error([' *** ERROR in FEM_PLAN:  An input structure with', ...
         ' sagittal femoral segmentations is required!']);
end
%
if (nargin<2)||isempty(iplt)
  iplt = false;
else
  iplt = true;
end
%
if (nargin<3)||isempty(ltxt)
  ltxt = [];
end
%
if (nargin<4)||isempty(fsav)
  iprt = false;
else
  iprt = true;
end
%
% Convert Radians to Degrees
%
rad2deg = 180/pi;
%
% Posterior and Medial/Lateral Cutoffs
%
% tmin = -145/rad2deg;    % Angle (theta) cutoff (-150, -145, or -140)
tmin = -140/rad2deg;    % Angle (theta) cutoff (-150, -145, or -140)
y0 = 0;                 % Y cutoff (-1, 0, or 1)
hpi = pi/2;
%
% Lateral-Medial Plane Normal Vector
%
lmp = [0 0 0];          % Lateral-medial point in the plane
lmv = [0 1 0];          % Lateral-medial plane normal vector
%
if iplt
  hf2 = figure;
  hold on;
  plt_datsl(f3,'b.-',0.5,[0 0 0.7]);
  axis equal;
  view(3);
  grid on;
  xlabel('X (mm)','FontSize',12,'FontWeight','bold');
  ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
  zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
  title({'Dividing Planes in Femur CS'; ltxt},'FontSize',16, ...
        'FontWeight','bold');
  view(40,20);
  axis equal;
  %
  axlim = axis;
  pxyz = [axlim(1) 0 axlim(5); axlim(2) 0 axlim(5); ...
          axlim(2) 0 axlim(6); axlim(1) 0 axlim(6); axlim(1) 0 axlim(5)];
  patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),pxyz(:,3),'FaceColor', ...
        [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
end
%
% Convert XZ to Polar Coordinates
%
xyzlbf = cell2mat(f3);
%
[thlf,rlf,zlf] = cart2pol(xyzlbf(:,1),xyzlbf(:,3),xyzlbf(:,2));
idc = thlf>hpi;
thlf(idc) = thlf(idc)-2*pi;
%
if iplt
  figure;
  plot3(thlf,zlf,rlf,'.');
  axlim = axis;
  view(2);
  hold on;
  plot([tmin tmin],axlim(3:4),'r-','LineWidth',1);
  xlabel('Theta (radians)','FontSize',12,'FontWeight','bold');
  ylabel('Z (mm)','FontSize',12,'FontWeight','bold');
  title({'Theta-Z CS'; ltxt},'FontSize',16,'FontWeight','bold');
else
  tzrmn = min([thlf,rlf,zlf]);
  tzrmx = max([thlf,rlf,zlf]);
  axlim = [tzrmn; tzrmx];
  axlim = axlim(:)';
end
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
%
ppc = mean(ppxyz(1:4,:));              % Center of plane
[nx,ny,nz] = plnorm(ppxyz(1:3,1),ppxyz(1:3,2),ppxyz(1:3,3));    % Normal
ppn = [nx,ny,nz];
%
if (iplt)
  figure(hf2);
  patch(ppxyz(:,1),ppxyz(:,2),ppxyz(:,3),ppxyz(:,3),'FaceColor', ...
        [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
end
%
% Get Maximum "Peaks" of Femur in X-Direction
%
nsl = size(f3,1);       % Number of slices
%
mx = zeros(nsl,1);
my = zeros(nsl,1);
%
for ks = 1:nsl
   [mx(ks),idx] = max(f3{ks}(:,1));
   my(ks) = f3{ks}(idx,2);
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
% Get Lateral and Medial Trochlear Planes
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
ppl = [0 0 0];          % Point in plane
%
zvm  = cross(xv,vm);
pnm = cross(vm,zvm);    % Medial plane normal
pnm = pnm./norm(pnm);
ppm = [0 0 0];          % Point in plane
%
% Plot Lateral and Medial Trochlear Planes
%
if iplt
%
  figure(hf2);
  orient landscape;
  axis equal;
  axlim = axis;
%
  xlmx = -pnl(2)*axlim(4)/pnl(1);
  xmmx = -pnm(2)*axlim(3)/pnm(1);
%
  pxyztll  = [0 0 axlim(6); 0 0 axlim(5); xlmx axlim(4) axlim(5); ...
              xlmx axlim(4) axlim(6); 0 0 axlim(6)];
%
  pxyztlm  = [0 0 axlim(6); xmmx axlim(3) axlim(6); xmmx ...
              axlim(3) axlim(5); 0 0 axlim(5); 0 0 axlim(6)];
%
  patch(pxyztll(:,1),pxyztll(:,2),pxyztll(:,3),pxyztll(:,3), ...
        'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
  patch(pxyztlm(:,1),pxyztlm(:,2),pxyztlm(:,3),pxyztlm(:,3), ...
        'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
%
  if iprt
    print('-dpsc2','-r600','-fillpage',fsav);
    view(-90,90);
    print('-dpsc2','-r600','-fillpage','-append',fsav);
  end
%
end
%
return