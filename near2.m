function [idx,xy1] = near2(xy0,xy)
%NEAR2     Finds the nearest point in a list nearest a reference point.
%
%          IDX = NEAR2(XY0,XY) finds the point in a list nearest a
%          reference point.  The list of points, XY, is an array of
%          X coordinates in the first column, Y coordinates in the
%          second column, and optional columns for additional
%          coordinates.  The reference point, XY0, is a row vector with
%          coordinates in the columns.  The index, IDX, into the list
%          to the point nearest the reference point is returned.
%
%          [IDX,XY1] = NEAR2(XY0,XY) returns the coordinates of the
%          nearest point in row vector XY1.
%
%          NOTES:  1.  Returns the first point nearest the reference
%                  point if the points are not unique.
%
%        24-Jan-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in NEAR2:  Two input arguments are required!');
end
%
% Calculate Squared Distances to Reference Point
%
npts = size(xy,1);      % Number of points in list
d = xy-repmat(xy0,npts,1);   % Differences
d = d.*d;               % Squared
d = sum(d,2);           % Sum of squared differences
[~,idx] = min(d);       % Find index to minimum
%
xy1 = xy(idx,:);        % Coordinates of minimum
%
return
