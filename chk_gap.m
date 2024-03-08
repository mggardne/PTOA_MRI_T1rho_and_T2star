function dat = chk_gap(dat)
%CHK_GAP   Checks distance between 2-D or 3-D points in a cellarray and
%          ensures the greatest distance is between the first and last
%          points.
%
%          DAT = CHK_GAP(DAT) given a cell array, DAT, containing two
%          (2) or three (3) columns matrices with coordinate point
%          data, returns a column cell array, DAT, with the points
%          within each cell reordered so the greatest distance between
%          points is between the first and last points.
%
%          NOTES:  1.  Every cell must contain at least two or more
%                  points.
%
%                  2.  M-file fix_gap.m must be in the current path or
%                  directory.
%
%          07-Mar-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Input
%
if (nargin<1)
  error(' *** ERROR in CHK_GAP:  Must have an input cell array!');
end
%
if ~iscell(dat)
  error(' *** ERROR in CHK_GAP:  Input data must be a cell array!');
end
%
% Check for Gaps in Point Coordinates
%
dat = cellfun(@fix_gap,dat(:),'UniformOutput',false);
%
return