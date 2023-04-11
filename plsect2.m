function [ip,t,idx] = plsect2(pp,pn,pwl,tol)
%PLSECT2 Finds the intersection of a plane and a piecewise linear (PWL)
%        line.
%
%        IP = PLSECT2(PP,PN,PWL) finds the first intersection of a
%        plane defined by a point on the plane, PP, and a vector normal
%        to the plane, PN, and a piecewise linear (PWL) line.  PWL is
%        defined by a series of 3-D points with the X, Y and Z
%        coordinates of the points in columns.  The X, Y, and Z
%        coordinates of the first intersection with the piecewise
%        linear line are returned in IP.  IP is an empty array if there
%        is no intersection.
%
%        [IP,T,IDX] = PLSECT2(PP,PN,PWL) returns the distance, T, along
%        the line with the first intersection.  T is the normalized
%        (0 to 1) distance along the line to the intersection.  The
%        index, IDX, is to the first point of the segment in the
%        piecewise linear line (PWL) that intersects the plane.  If
%        there is no intersection, T and IDX are empty arrays.
%        
%        IP = PLSECT2(PP,PN,PWL,TOL) checks to make make sure the dot
%        product of the line directions and plane normal is not within
%        tolerance TOL of zero indicating the line is (or close to)
%        parallel to the plane.  Default tolerance is 1e-8.
%        
%        NOTES:  1.  The M-file plsect.m must be in the current path or
%                directory.
%
%       02-Dec-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** Error in PLSECT2:  Three input arguments are required!');
end
%
if (nargin<4)||isempty(tol)
  tol = 1e-8;
end
%
% Get Column Vectors
%
pp = pp(:);
pn = pn(:);
%
% Check Inputs
%
[npts,nc] = size(pwl);
if size(pp,1)~=3||size(pn,1)~=3||nc~=3
  error([' *** Error in PLSECT2:  Coordinate dimensions must equal', ...
         ' three (3)!']);
end
if npts<2
  error([' *** Error in PLSECT2:  Piecewise linear line must have', ...
         ' at least two (2) points!']);
end
%
% Initial Values
%
ip = [];                % Empty if no intersection
t = [];                 % Empty if no intersection
idx = [];               % Empty if no intersection
%
nl = npts-1;            % Number of lines in piecewise linear (PWL) line
%
% Loop through Piecewise Linear Line
%
for k = 1:nl
%
   l = k+1;             % Next point
%
   lp = pwl(k,:)';                     % Point on line
   lv = pwl(l,:)'-lp;                  % Direction of line
%
% Check for Intersection
%
   [ipp,tp] = plsect(pp,pn,lp,lv,tol);
%
   if ~isempty(ipp)
     ip = ipp;
     t = tp;
     idx = k;
     break;
   end
%
end
%
return