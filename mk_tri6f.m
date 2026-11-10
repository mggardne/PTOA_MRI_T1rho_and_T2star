function [tri,nt,slx] = mk_tri5f(dat,ang,tol,iplt)
%MK_TRI5F Makes a triangular mesh by using the ordered slice data from
%         the digitized MRI femur data.  Checks the anterior-posterior
%         coordinates of adjacent slices to ensure that they overlap.
%         If there is no overlap, the first two points of one slice and
%         the last two points of the next slice are connected by a pair
%         of triangles. 
%
%         [TRI,NT] = MK_TRI5F(DAT) given a cell array containing three
%         (3) columns matrices with slice coordinate point data, DAT,
%         returns the three (3) column triangle connectivity matrix,
%         TRI.  The number of returned triangles, NT, may also be
%         returned.
%
%         TRI = MK_TRI5F(DAT,ANG) The data is converted to cylindrical
%         coordinates with the Y-axis forming the Z-axis of the cylinder
%         (X-Z plane transformed to polar coordinates). The mesh is only
%         defined between the angles of the smaller of the two slices
%         plus an additional angle of ANG degrees.  By default, ANG is
%         three (3) degrees.
%
%         NOTES:  1.  Each slice coordinate data matrix must correspond
%                 to one index into the cell array DAT.
%
%                 2.  The coordinates should be ordered in the same
%                 direction in every slice.  The dot product of the
%                 directions of adjacent slices are used to check the
%                 ordering direction and the ordering direction is
%                 reversed if the dot product is negative.
%
%                 3.  The arclength along each slice is used to
%                 determine the triangulation.  The data is converted
%                 to cylindrical coordinates with the Y-axis forming
%                 the Z-axis of the cylinder (X-Z plane transformed to
%                 polar coordinates).
% 
%                 4.  The anterior direction is assumed to be the +X
%                 direction (positive first coordinate).
% 
%                 5.  The M-files mk_tri4f.m, plane_fit.m, sl_info.m,
%                 and xzpl2pol.m must be in the current path or
%                 directory.
%
%         21-Apr-2025 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  iplt = false;
end
%
if iplt
  figure;
  orient landscape;
end
%
if (nargin<3)||isempty(tol)
  tol = 0.1;            % Tolerance on dot product (projection of directions)
end
%
if (nargin<2)||isempty(ang)
  ang = 3;
end
%
if (nargin<1)
  error(' *** ERROR in mk_tri5f:  No input data!');
end
%
% Check that Adjacent Slices Overlap
%
dat = dat(:);
nsl = size(dat,1);
nslm1 = nsl-1;
%
dx = zeros(nslm1,1);
%
for k = 1:nsl-1
   dx(k) = dat{k+1}(end,1)-dat{k}(1,1);
end
%
% Check for No Overlap in the Anterior-Posterior Coordinates
%
idx = find(dx>0);
%
% If All Adjacent Slices Overlap, Use mk_tri4f Directly
%
if isempty(idx)
  [tri,nt,slx] = mk_tri4f(dat,ang,tol,iplt);
else
%
% Loop Through Slices Which Do Not Overlap
%
  nx = size(idx,1);
  idx = [idx idx+1];    % Adjacent slices with no overlap
  idx = [0 1; idx; 42 0];
  [~,~,isl] = sl_info(dat);            % Index to slice coordinates
%
  tri = cell(nx,1);
%
  for k = 1:nx
%
     idx1 = idx(k,2);
     idx2 = idx(k+1,1);
     tri{k} = mk_tri4f(dat(idx1:idx2),ang,tol,iplt);
     if k>1
       tri{k} = tri{k}+isl(idx(k+1,2));
     end
     pids = [isl(idx2)+1 isl(idx2)+2 isl(idx2+2)-1 isl(idx2+2)];
     tri{k} = [tri{k}; pids(3) pids(1) pids(4); pids(1:2) pids(4)];
%  
  end
%
  k = k+1;
  idx1 = idx(k,2);
  idx2 = idx(k+1,1);
  trit = mk_tri4f(dat(idx1:idx2),ang,tol,iplt);
  trit = trit+isl(idx1);
  tri = [cell2mat(tri); trit];
%
end
%
return