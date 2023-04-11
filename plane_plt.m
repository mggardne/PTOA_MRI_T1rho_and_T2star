function [xp,yp,zp] = plane_plt(pxyz,nv,mnmx,na);
%PLANE_PLT X, Y and Z mesh coordinate data for a plane.
%
%         [XP,YP,ZP] = PLANE_PLT(PXYZ,NV,MNMX) given a vector with a
%         point (X, Y and Z) on a plane, PXYZ, the normal vector, NV,
%         and a matrix with the minimum extent of the plane in the
%         first row and the maximum extent of the plane in the second
%         row for all three coordinates (X, Y and Z) in the columns,
%         returns the coordinate matrices, XP, YP and ZP, defining the
%         surface of a plane.
%
%         [XP,YP,ZP] = PLANE_PLT(PXYZ,NV,MNMX,NA) NA is a two (2)
%         element vector of the number of points in each direction and
%         returns the coordinates of NA points in each direction of the
%         plane.  Generally, NA >= 5 with the default of NA = 10.  Use
%         NA = 2 to just get the outer boundary of the plane.
%
%         NOTES:  1.  See PLANE_FIT.M for fitting a plane.
%
%         24-June-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error([' *** ERROR in PLANE_PLT:  At least three (3) inputs', ...
         ' are required!']);
end
%
if (nargin<4)
  na = [10; 10];
end
%
na = na(:);
if size(na,1)==1
  na = [na; na];        % Allow different mesh density in each direction
end
%
% Find Maximum Direction of Normal and Index for Other Two Coordinates
%
pxyz = pxyz(:)';        % Row vector
nv = nv(:);             % Column vector
%
[~,idmx] = max(abs(nv));
idc = 1:3;
idc(idmx) = 0;
idc = find(idc);
%
mnmx = mnmx(:,idc);
%
% Get Coordinates
%
if idmx==1
  [yp,zp] = meshgrid(linspace(mnmx(1,1),mnmx(2,1),na(1)), ...
            linspace(mnmx(1,2),mnmx(2,2),na(2)));
  xp = 1/nv(idmx)*(pxyz*nv-yp*nv(2)-zp*nv(3));
elseif idmx==2
  [xp,zp] = meshgrid(linspace(mnmx(1,1),mnmx(2,1),na(1)), ...
            linspace(mnmx(1,2),mnmx(2,2),na(2)));
  yp = 1/nv(idmx)*(pxyz*nv-xp*nv(1)-zp*nv(3));
else
  [xp,yp] = meshgrid(linspace(mnmx(1,1),mnmx(2,1),na(1)), ...
            linspace(mnmx(1,2),mnmx(2,2),na(2)));
  zp = 1/nv(idmx)*(pxyz*nv-xp*nv(1)-yp*nv(2));
end
%
return