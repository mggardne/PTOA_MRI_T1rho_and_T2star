function [xyz,idx,npts] = fix_gap(xyz)
%FIX_GAP   Checks distance between 2-D or 3-D points and ensures the
%          greatest distance is between the first and last points.
%
%          XYZ = FIX_GAP(XYZ) given a two (2) or three (3) columns
%          matrix with coordinate point data, XYZ, returns the matrix,
%          XYZ, with the points reordered so the greatest distance
%          between points is between the first and last points.
%
%          [XYZ,IDX] = FIX_GAP(XYZ) The index, IDX, is returned such
%          that the returned XYZ = the input XYZ(IDX,:).
%
%          [XYZ,IDX,NPTS] = FIX_GAP(XYZ) The number of points, NPTS,
%          may also be returned.
%
%          NOTES:  None.
%
%          31-May-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Input
%
if (nargin<1)
  error(' *** ERROR in FIX_GAP:  Must have input points coordinates!');
end
%
[npts,ncol] = size(xyz);
%
if npts<2
  error(' *** ERROR in FIX_GAP:  Data must have at least two points!');
end
%
if ncol~=2&&ncol~=3
  error(' *** ERROR in FIX_GAP:  Data must have two or three columns!');
end
%
% Find Distances and Maximum Distance
%
d = diff([xyz; xyz(1,:)]);
d = sum(d.*d,2);        % Distances squared
%
[~,idmx] = max(d);      % Location of maximum distance
%
% Check Maximum Distance Location and Reorder Points if Necessary
%
idx = (1:npts)';
if idmx~=npts
  idx = [idmx+1:npts 1:idmx]';
  xyz = xyz(idx,:);
end
%
return