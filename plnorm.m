function [nx,ny,nz,xc,yc,zc] = plnorm(x,y,z)
%PLNORM Computes the normal and centroid of a 3-D triangle.
%       [Nx,Ny,Nz,Xc,Yc,Zc] = PLNORM(X,Y,Z) returns the normal and
%       centroid of a triangle defined by three (3) points with
%       coordinates X,Y and Z.  X, Y and Z must be vectors of length
%       three (3).  The coordinates should be listed in clockwise
%       direction when seen from above.
%
%       28-Mar-96

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error('PLNORM requires three input arguments.');
end
%
% Check that Inputs are Vectors
%
x = x(:);
y = y(:);
z = z(:);
%
[n l] = size(x);
[n2 l2] = size(y);
[n3 l3] = size(z);
%
if ((l~=1)|(l2~=1)|(l3~=1))
  error('PLNORM only works with vectors.')
end
%
% Check that the Inputs have Three Rows
%
if ((n~=3)|(n2~=3)|(n3~=3))
  error('X, Y and Z must be of length three (3).')
end
%
% Normal of Triangle
%
a = [(x(1)-x(2)); (y(1)-y(2)); (z(1)-z(2))];
b = [(x(3)-x(2)); (y(3)-y(2)); (z(3)-z(2))];
n = cross(a,b);
n = n/norm(n);
%
nx = n(1);
ny = n(2);
nz = n(3);
%
% Centroid of Triangle
%
xc = mean(x);
yc = mean(y);
zc = mean(z);
%
return