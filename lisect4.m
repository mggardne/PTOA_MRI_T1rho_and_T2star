function [ip,t,idx,ierr] = lisect4(p1,v1,pwl,tol,clim)
%LISECT4  Finds the first intersection of a closed line with a piecewise
%         linear line.
%
%         IP = LISECT4(P1,V1,PWL) finds the first intersection of a
%         closed line defined by a point (P1) and a vector (V1) with a
%         piecewise linear (PWL) line.  PWL is defined by a series of
%         3-D points with the X, Y and Z coordinates of the points in
%         columns. The X, Y and Z coordinates of the first intersection
%         with the piecewise linear line are returned in IP.  IP is
%         empty if there is no intersection.
%
%         [IP,T,IDX] = LISECT4(P1,V1,PWL) returns the distance along the
%         lines to the intersection.  T(1) is the normalized (0 to 1)
%         distance along the first line to the intersection.  T(2) is
%         the distance along the segment of the piecewise linear line
%         (PWL) with the intersection.  The index (IDX) is to the first
%         point of the segment in the piecewise linear line (PWL) that
%         intersects the input line.  If there is no intersection, T
%         and IDX are empty arrays.
%
%         [IP,T,IDX,IERR] = LISECT4(P1,V1,PWL,TOL) sets IERR to true if
%         no intersection is found within tolerance (TOL).  Default
%         tolerance is 1e-8.
%
%         [IP,T,IDX,IERR] = LISECT4(P1,V1,PWL,TOL,CLIM) sets IERR to
%         true if the condition number of the matrix is greater than a
%         limit, CLIM.  This is usually due to parallel or nearly
%         parallel lines.  The default condition number limit is 1e+8.
%
%         NOTES:  1.  The lines are assumed to be not parallel.
%
%                 2.  The input point and vector must be of length
%                 three (3).  The piecewise linear line must have three
%                 (3) columns.
%
%                 3.  See lsect.m and lsect2.m for two-dimensional
%                 line (2-D) intersections.
%
%                 4.  See lisect.m and lisect2.m for three-dimensional
%                 line (3-D) intersections with infinite lines.  See
%                 lisect3.m for the intersection between two closed
%                 lines.
%
%         23-May-2014
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error([' *** ERROR in  LISECT4:  LISECT4 requires three input', ...
         ' arguments.']);
end
%
if (nargin==3)|(tol<=0)
  tol = 1e-8;
end
%
if (nargin<4)||(isempty(tol))||(tol<=0)
  tol = 1e-8;
end
%
if (nargin<5)||(isempty(clim))||(clim<=1)
  clim = 1e+8;
end
%
% Check Vectors and Piecewise Linear (PWL) Line
%
p1 = p1(:);
v1 = v1(:);
[n1 l1] = size(p1);
[n2 l2] = size(v1);
[n3 l3] = size(pwl);
%
if ((n1~=3)|(l1~=1)|(n2~=3)|(l2~=1)|(n3<2)|(l3~=3))
  error(' *** ERROR in LSECT4:  Error in input arguments.')
end
%
% Initialize Arrays
%
nl = n3-1;              % Number of lines in piecewise linear (PWL) line
%
ip = [];
t = [];
idx = [];
ierr = false;
%
% Loop through Piecewise Linear Line
%
for k = 1:nl
%
   l = k+1;             % Next point
%
   p2 = pwl(k,:)';                     % Point on line
   v2 = pwl(l,:)'-p2;                  % Direction of line
%
% Intersection Equations
%
   A = [v1(1:2) -v2(1:2)];             % Solve for intersection
   b = p2(1:2)-p1(1:2);
%
   if (rank(A)~=2)
     A = [v1(2:3) -v2(2:3)];
     b = p2(2:3)-p1(2:3);
   end
%
   if (rank(A)~=2)
     A = [v1([1 3]) -v2([1 3])];
     b = p2([1 3])-p1([1 3]);
   end
%
% Check Condition Number
%
   cn = cond(A);
   if cn>clim
     warning([' *** WARNING in LISECT4:  Lines are parallel' ...
              ' within the condition number limit.']);
     ierr = true;
     return
   end
%
% Solve for Intersection
%
   tv = A\b;
%
   xyz1 = tv(1)*v1+p1;                 % Coordinates of intersection
   xyz2 = tv(2)*v2+p2;
%
   chk = norm(xyz1-xyz2);              % Check intersection
   if (chk>tol)
     warning([' *** WARNING in LISECT4:  Lines do not intersect', ...
              ' within tolerance.']);
     ierr = true;
     return
   end
%
% Check Within Closed Lines
%
   if all(tv>=0)&&all(tv<=1)           % Incorporate tolerance?
     ip = xyz1';
     t = tv;
     idx = k;
     return;
   end
%
end
%
return