function datt = trim_seg(dat,dist,tol,iplt)
%TRIM_SEG  Trims a second segmentation based on a first segmentation.
%          Both segmentations must lay in the same plane.
%
%          DATT = TRIM_SEG(DAT) given a two-element cell array
%          containing two (2) segmentations with three-dimesional (3D)
%          coordinate point data, DAT, returns the two (2)
%          segmentations with the second segmentation trimmed so its
%          length is close to the first segmentation in a two-element
%          cell array DATT.
%
%          DATT = TRIM_SEG(DAT,DIST) given a distance, DIST, ensures
%          truncated endpoints on the second segmentation line are
%          within two DIST distances of the first segmentation line
%          endpoints.
%
%          NOTES:  1.  Segmentation coordinates for both lines must lay
%                  within the same plane.
%
%                  2.  Each coordinate data matrix must correspond to
%                  one index into the cell array DAT.
%
%                  3.  The coordinates should be ordered in the same
%                  direction in each line.  The dot product of the
%                  directions of the adjacent lines are used to check
%                  the ordering direction and the ordering direction is
%                  reversed if the dot product is negative.
%
%                  4.  The second segmentation line may be truncated if
%                  it extends pass the first segmentation line.
%
%                  5.  M-files lsect2a.m, plane_fit, and near2.m must
%                  be in the current directory or path.
%
%          25-Jan-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in trim_seg:  No input data!');
end
%
if nargin<2||isempty(dist)
  dist = Inf;            % No distance checking
end
%
diste = 2*dist;         % Distance from cartilage end points to line of bone points
diste = diste*diste;    % Squared distance for comparison to squared lengths
%
if nargin<3||isempty(tol)
  tol = 0.01;
end
%
if nargin<4||isempty(iplt)
  iplt = false;
end
%
dat = dat(:);
nslice = size(dat,1);
%
if nslice~=2
  error([' *** ERROR in trim_seg:  Input cell array must have two', ...
         ' elements containing coordinates for two lines!']);
end
%
% Get 3D Data and Convert to 2D Data
%
xyz1 = dat{1};          % Cartilage coordinates
npts1 = size(xyz1,1);
xyz2 = dat{2};          % Bone coordinates
npts2 = size(xyz2,1);
%
xyz = [xyz1; xyz2];
npts = npts1+npts2;
[~,~,~,r] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));    % Rotation matrix
%
xyz0 = mean(xyz);       % Origin
xy = xyz-repmat(xyz0,npts,1);          % Center
xy = xy*r;              % Rotate into plane
xy = xy(:,1:2);
%
xy1 = xy(1:npts1,:);
xy2 = xy(npts1+1:npts,:);
%
% Get First (Cartilage) Line and Slopes at the Ends of the Line
%
vec1 = xy1(npts1,:)-xy1(1,:);          % Direction of line
vec1 = vec1./norm(vec1);
%
dd = diff(xy1);
de = dd([1,npts1-1],:);                % Directions at ends of first (top) line
%
% Equations of Lines at the Ends of the First (Cartilage) Line
%
mp = -de(:,1)./de(:,2); % 90 degrees
ids = abs(de(:,2))<1e-8;
mp(ids) = sign(mp(ids))*1e+4;          % Large but not infinite slope
% mp(2,1) = (de(2,1)+de(2,2))./(de(2,1)-de(2,2));  % 45 degrees
% mp(1,1) = (de(1,2)-de(1,1))./(de(1,1)+de(1,2));  % 45 degrees
% ratio = tan(pi/3);      % Opposite to adjacent for 60 degrees (pi/3 radians)
% mp(2,1) = (de(2,2)+ratio*de(2,1))./(de(2,1)-ratio*de(2,2));     % 60 degrees
% mp(1,1) = (ratio*de(1,2)-de(1,1))./(de(1,1)+ratio*de(1,2));     % 60 degrees
xp = xy1([1; npts1],1);
yp = xy1([1; npts1],2);
bp = yp-mp.*xp;
%
% Get Direction of Second (Bone) Line
%
vec2 = xy2(npts2,:)-xy2(1,:);
vec2 = vec2./norm(vec2);
%
% Check for Slices with a Reverse Digitization
%
dotp = vec1*vec2';
%
if dotp<tol
  xy2 = flipud(xy2);
  vec = xy2(npts2,:)-xy2(1,:);
  vec = vec./norm(vec);
  dotp2 = vec1*vec';
  if dotp2<dotp         % Revert back to original ordering
    warning([' *** WARNING in trim_seg:  Ordering of points', ...
             ' in the slices may not be in the same direction!']);
    xy2 = flipud(xy2);
  end
end
%
% Cut Off Extra Points on Second (Bone) Line
% Fails if Second (Bone) is Not Longer Than First (Cartilage) Line
%
[xyi,id1] = lsect2a(mp(1),bp(1),xy2);
ni = size(id1,1);
if ni>1
  idp = near2([xp(1) yp(1)],xyi);      % Get nearest point
  id1 = id1(idp);
  xyi = xyi(idp,:);
end
d = xyi-xy1(1,:);
d = d*d';
if d>diste              % Bone line within 2 DIST of cartilage end points
  id1 = 2;
end
if isempty(id1)
  id1 = 1;
end
%
[xyi,id2] = lsect2a(mp(2),bp(2),xy2);
ni = size(id2,1);
if ni>1
  idp = near2([xp(2) yp(2)],xyi);      % Get nearest point
  id2 = id2(idp);
  xyi = xyi(idp,:);
end
d = xyi-xy1(npts1,:);
d = d*d';
if d>diste              % Bone line within 2 DIST of cartilage end points
  id2 = npts2-1;
end
if isempty(id2)
  id2 = npts2-1;
end
%
idc = id1:id2+1;
xy2 = xy2(idc,:);
npts2_new = size(xy2,1);
%
% Convert Back to 3D
%
xyz2 = [xy2 zeros(npts2_new,1)]*r'+repmat(xyz0,npts2_new,1);
%
datt = dat;
datt{2} = xyz2;
%
if iplt
  lincolr = ['b.-'; 'g.-'];
  figure;
  hold on;
  orient landscape;
  for k = 1:2
     plot3(dat{k}(:,1),dat{k}(:,2),dat{k}(:,3),lincolr(k,:), ...
           'LineWidth',2,'MarkerSize',8);
  end
  plot3(xyz2(:,1),xyz2(:,2),xyz2(:,3),'r.-','LineWidth',1, ...
        'MarkerSize',7);
  view(25,20);
  axis equal;
end
%
return