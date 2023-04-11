function [ip,t,il,ierr] = plsect(pp,pn,lp,lv,tol)
%PLSECT  Finds the intersection of a plane and a line.
%
%        [IP,T,IL,IERR] = PLSECT(PP,PN,LP,LV) finds the intersection of
%        a plane defined by a point on the plane, PP, and a vector
%        normal to the plane, PN, and a line defined by a
%        point LP and a direction vector LV.  The lengths of all the
%        inputs must be three (3) for the X, Y, and Z components of the
%        coordinates or vectors.  IP is the coordinates of the
%        intersection, T is the parametric distance along the line to
%        the intersection.  IL is set to true (1) if an intersection
%        was found.  IERR is set to true (1) if the line is (or close
%        to) parallel to the plane.  If there is no intersection, IP
%        and T are empty arrays.
%        
%        IP = PLSECT(PP,PN,LP,LV,TOL) checks to make make sure the dot
%        product of the line direction and plane normal is not within
%        tolerance TOL of zero indicating the line is (or close to)
%        parallel to the plane.  Default tolerance is 1e-8.
%        
%        NOTES:  None.
%
%       02-Dec-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  error(' *** Error in PLSECT:  Four input arguments are required!');
end
%
if (nargin<5)||isempty(tol)
  tol = 1e-8;
end
%
% Get Column Vectors
%
pp = pp(:);
pn = pn(:);
lp = lp(:);
lv = lv(:);
%
% Initial Values
%
ip = [];                % Empty if no intersection
t = [];                 % Empty if no intersection
il = false;
ierr = false;
%
% Check Dot Product of the Line Direction and Plane Normal
%
if abs(lv'*pn)<tol
  ierr = true;
  return;
end
%
% Calculate Intersection Point
% Solution of substituting the parametric line equation into the plane
% equation.
%
t = pn'*(pp-lp)./(pn'*lv);             % Parametric distance along line
%
% Check Intersection is Within the Line Endpoints (0<=t<=1)
%
if (0<=t)&&(t<=1)
  il = true;
  ip = lp+t*lv;
  ip = ip';             % Coordinates in columns
end
%
return