function [xyzp,t] = pt2line(lpt0,lvec,pxyz);
%PT2LINE  Calculates the coordinates of a point on a line that is
%         closest to a point not on the line.  This point also forms a
%         perpendicular to the line that passes through the point not
%         on the line.
%
%         XYZP = PT2LINE(LPT0,LVEC,PXYZ) given a line defined by an
%         initial point on the line in a row vector, LPT0, with X, Y
%         and Z coordinates in columns, a direction vector for the line
%         in a row vector, LVEC, with  X, Y and Z components and the
%         coordinates of the point not on the line in a row vector,
%         PXYZ, with X, Y and and Z coordinates in columns, calculates
%         the coordinates of a point on the line that is closest to
%         the point not on the line.  This point also forms a
%         perpendicular to the line that passes through the point not
%         on the line.  The X, Y and Z coordinates of the point on the
%         line is returned in the row vector, XYZP.
%
%         [XYZP,T] = PT2LINE(LPT0,LVEC,PXYZ) returns the parametric
%         distance, T, along the line to the point on the line that is
%         closest to the point not on the line.
%
%         NOTES:  1.  The line is assumed to have infinite length.
%
%                 2.  All points and the vector must have three columns.
%
%         20-July-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error([' *** ERROR in PT2LINE:  At least three (3) inputs', ...
         ' are required!']);
end
%
% Check Inputs
%
[nrow1,ncol1] = size(lpt0);
[nrow2,ncol2] = size(lvec);
[nrow3,ncol3] = size(pxyz);
if nrow1~=1|nrow2~=1|nrow3~=1|ncol1~=3|ncol2~=3|ncol3~=3
  error([' *** ERROR in PT2LINE:  All input points and the vector', ...
         ' must have three columns!']);
end
%
% Coordinates of the Perpendicular (Closest) Point on the Line
%
dp = pxyz-lpt0;                        % Vector between point and beginning of line
t = (lvec*dp')/(lvec*lvec');           % Parametric distance along line to perpendicular point
xyzp = t*lvec+lpt0;                    % Coordinates of the perpendicular point on the line
%
return