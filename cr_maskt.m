function mask = cr_maskt(roic,npx,dist,tol,iplt)
%CR_MASKT  Creates an image mask based on the two-dimensional (2-D)
%          coordinates of two lines.  The region of interest (ROI) is
%          defined as the space between the first and second line.  The
%          second line is assumed to be slightly longer than the first
%          line.
%
%          For creating masks for knee tibia cartilage.  The first line
%          is usually cartilage and the second longer line is the
%          underlying bone.
%
%          MASK = CR_MASKT(ROIC,NPX) given the two-dimensional
%          coordinates of two lines defining a region of interest in a
%          two element cell array (one element for each line), ROIC,
%          and the size of the image for the mask in array, NPX,
%          creates a logical array mask for the region of interest,
%          MASK.  Note:  If NPX has only one element, the image is
%          assumed to be square (symmetrical) (NPX by NPX).
%
%          NOTES:  1.  M-files in_tri2d.m, lsect2a.m, near2.m, and
%                  mk_tri_2d.m must be in the current directory or
%                  path.
%
%                  2.  Similar to cr_maskf.m except the slopes of the
%                  end points are in the opposite directions.
%
%          04-Nov-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in CR_MASKT:  Two inputs are required!');
end
%
roic = roic(:);
[nr,nc] = size(roic);
if nr~=2&&nc~=1
  error([' *** ERROR in CR_MASKT:  Input cell array must have'
         ' two elements!']);
end
%
ndim = size(npx(:),1);
if ndim>2||ndim<1
  error(' *** ERROR in CR_MASKT:  Incorrect image dimensions!');
end
%
if ndim==1
  npx = [npx; npx];
end
%
if nargin<3||isempty(dist)
  dist = Inf;           % No distance checking
end
%
if nargin<4||isempty(tol)
  tol = 0.1;
end
%
if nargin<5
  iplt = false;
end
%
% Create Triangle Meshes for the Region of Interest
%
[tri,xy] = mk_tri_2d(roic,dist,tol,iplt);
%
% Find Pixels within the Region of Interest
%
mask = false(npx(1)*npx(2),1);
%
minr = floor(min(xy));
if any(minr(:)==0)      % Trap for zeros
  minr(minr(:)==0) = 1;
end
maxr = ceil(max(xy));
idx = minr(:,1):maxr(:,1);
idy = minr(:,2):maxr(:,2);
[xg,yg] = meshgrid(idx,idy);
xym = [xg(:) yg(:)];
in_roi = in_tri2d(tri,xy,xym);
%
idr = sub2ind([npx(1) npx(2)],xym(:,2),xym(:,1));
idr = idr(in_roi);
%
mask(idr) = true;
%
return