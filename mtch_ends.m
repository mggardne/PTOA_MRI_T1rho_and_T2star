function [varargout] = mtch_ends(tr,tol);
%MTCH_ENDS Checks rows of a cell array of coordinates to see if the
%          first or last points from the first column have the
%          same coordinates as the first or last points in the other
%          columns.
%
%          To control the clipping of the bone segmentation lines when
%          the bone and cartilage are broken into regions of interest.
%
%          [MTCH1,MTCH2, ...] = MTCH_ENDS(TR,TOL) given a cell array
%          with NROW rows and NCOL columns with sets of coordinates, TR,
%          checks the rows of TR to see if the first or last points from
%          the first column have the same coordinates within a distance
%          tolerance, TOL, as the first or last points in the other
%          columns.  Four (4) by NCOL-1 logical arrays are returned for
%          each row.
%
%          NOTES:  1.  Default tolerance is 0.1.
%
%          07-Feb-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Check Inputs
%
if (nargin<1)||~iscell(tr)
  error(' *** ERROR in MTCH_ENDS:  An input cell array is required!');
end
%
if (nargin<2)
  tol = 0.1;
end
%
tol = tol*tol;          % Square for comparison to the squared distances
%
[nr,ncol] = size(tr);
%
if ncol<2
  error([' *** ERROR in MTCH_ENDS:  At least two sets of', ... 
         ' coordinates in the columns of the cell array are'
         ' required!']);
end
%
% Loop through the Rows and Columns of the Cell Array
%
nc = ncol-1;            % Number of comparisons
%
npts = cellfun('length',tr);           % Number of points in each cell
%
id14 = [true; false; false; true];     % Index to rows of logical arrays
id23 = [false; true; true; false];     % Index to rows of logical arrays
%
varargout = cell(1,nr);                % Output logical arrays
%
for k = 1:nr
%
   npt = npts(k,:);     % Number of points in the columns of this row
   mtch = false(4,nc);
   pts0_14 = tr{k,1}([1; npt(1)],:);
   pts0_23 = pts0_14([2; 1],:);
%
   for l = 1:nc
%
      n = l+1;
      d = tr{k,n}([1; npt(n)],:);
%
      dd = d-pts0_14;
      dd = sum(dd.*dd,2);
      mtch(id14,l) = dd<tol;
%
      dd = d-pts0_23;
      dd = sum(dd.*dd,2);
      mtch(id23,l) = dd<tol;
%
   end
%
   varargout{k} = mtch;
%
end
%
return
      

