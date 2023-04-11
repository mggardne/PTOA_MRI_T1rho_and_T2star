function hc = cyl_plt(r,xyz1,xyz2,n,fc);
%CYL_PLT  Plots a cylinder given a radius and the centers of the two
%         ends of a cylinder.
%
%         HC = CYL_PLT(R,XYZ1,XYZ2) given the radius of the cylinder R,
%         a row vector of the coordinates of the center of one end of
%         the cylinder XYZ1, and a row vector of the coordinates of the
%         center of the other end of the cylinder XYZ2, plots the
%         cylinder as a lighted surface using SURFL and returns the
%         graphics handle HC.
%
%         HC = CYL_PLT(R,XYZ1,XYZ2,N) sets the number of circumferential
%         facets to N.  The default is 72.
%
%         HC = CYL_PLT(R,XYZ1,XYZ2,N,FC) sets the face color of the
%         cylinder to FC.  See "set(hc,'FaceColor')" for valid color
%         options.  By default there is no face color (only a wire
%         frame cylinder).
%
%         NOTES:  1.  For use with cyl_fit.m to plot the resulting
%                 cylinder.
%
%         17-Jul-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in CYL_PLT:  Three (3) inputs are required!');
end
%
if nargin<4||isempty(n)
  n = 72;
end
%
if nargin<5||isempty(fc)
  fc = 'none';
end
%
% Check Inputs
%
[nrow1,ncol1] = size(xyz1);
[nrow2,ncol2] = size(xyz2);
if nrow1~=1|nrow2~=1|ncol1~=3|ncol2~=3
  error([' *** ERROR in CYL_PLT:  All input points must have', ...
         ' three columns!']);
end
%
% Get Cylinder Parameters
%
zc = xyz2-xyz1;
l = zc*zc';
l = sqrt(l);                           % Length of the cylinder
zc = zc./l;                            % Unit direction of cylinder axis
m = (xyz2+xyz1)/2;                     % Center of the cylinder
%
% Get Cylinder Surface
%
[x,y,z] = cylinder([r r r],n);
z = l*z-l/2;                           % Scale cylinder length
dim = size(x);
%
% Get Rotation Matrix and Translation Vector
%
yc = zeros(1,3);
[val,idm] = min(abs(zc));
yc(idm) = 1;                           % Arbitrary Y axis
xc = cross(yc,zc);
yc = cross(zc,xc);
rot = [xc; yc; zc];     % Local to global rotation matrix
%
a = [x(:)'; y(:)'; z(:)'];
%
a = rot'*a;                            % Rotate cylinder
a = a+repmat(m',1,prod(dim));          % Translate cylinder
%
x = reshape(a(1,:),dim);
y = reshape(a(2,:),dim);
z = reshape(a(3,:),dim);
%
% Plot Cylinder
%
hc = surfl(x,y,z);
set(hc,'FaceColor',fc);
set(hc,'LineWidth',0.5);
%
return