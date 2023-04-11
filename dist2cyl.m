function dist = dist2cyl(x,xyzp);
%DIST2CYL Finds the distances from a set of points and the surface of a
%         cylinder.
%
%         DIST = DIST2CYL(X,XYZP) given the seven (7) parameter vector
%         X that describes the cylinder and a matrix of points XYZP
%         with the X, Y and and Z coordinates in columns, returns the
%         distances between the points and the surface of a cylinder.
%
%         The seven (7) parameters in X are:
%         1.  radius of the cylinder
%         2.  X coordinate of first point on the axis of the cylinder
%         3.  Y coordinate of first point on the axis of the cylinder
%         4.  Z coordinate of first point on the axis of the cylinder
%         5.  X coordinate of second point on the axis of the cylinder
%         6.  Y coordinate of second point on the axis of the cylinder
%         7.  Z coordinate of second point on the axis of the cylinder
%
%         NOTES:  1.  For use with cyl_fit.m.  See cyl_fit.m for more
%                 details.
%
%                 2.  The M-file pts2lin.m must be in the current path
%                 or directory.
%
%         17-Jul-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in DIST2CYL:  Two (2) inputs are required!');
end
%
% Check Inputs
%
x = x(:);
nparam = size(x,1);
if nparam~=7
  error([' *** ERROR in DIST2CYL:  Seven parameters are required', ...
         ' to define the cylinder!']);
end
%
[np,ncol] = size(xyzp);
if ncol~=3
  error([' *** ERROR in DIST2CYL:  Points coordinate matrix', ...
         ' must have three columns!']);
end
%
% Get Cylinder Parameters
%
r = x(1);
lpt0 = x(2:4)';
lvec = x(5:7)'-lpt0;
lvec = lvec./sqrt(lvec*lvec');         % Unit vector
%
% Get Distances to Cylinder Surface
%
xyzc = pts2lin(lpt0,lvec,xyzp);        % Coordinates of points on the axis of the cylinder
%
dist = xyzc-xyzp;
dist = sqrt(sum(dist.*dist,2));        % Distances to axis of the cylinder
dist = dist-r;                         % Distances to surface of the cylinder
%
return