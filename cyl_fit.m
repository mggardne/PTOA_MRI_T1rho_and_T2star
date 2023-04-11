function [r,xyz1,xyz2,sse,res,exitflag] = cyl_fit(ri,xyz1i,xyz2i, ...
                                                  xyzp,tol)
%CYL_FIT  Fits a cylinder to a set of points using nonlinear least
%         squares.
%
%         [R,XYZ1,XYZ2] = CYL_FIT(RI,XYZ1I,XYZ2I,XYZP) given an initial
%         radius RI, a row vector of the coordinates of an initial
%         point on the cylinder axis XYZ1I, an initial second point on
%         the cylinder axis XYZ2I, and a matrix of points XYZP with the
%         X, Y, and and Z coordinates in columns, returns the nonlinear
%         least squares fit of a cylinder with radius R, a point on the
%         axis of the cylinder XYZ1, and a second point on the axis of
%         the cylinder XYZ2.
%
%         [R,XYZ1,XYZ2] = CYL_FIT(RI,XYZ1I,XYZ2I,XYZP,TOL) set the
%         parameter and function tolerances in the optimization to TOL.
%
%         [R,XYZ1,XYZ2,SSE,RES] = CYL_FIT(RI,XYZ1I,XYZ2I,XYZP) returns
%         the sum of squared errors SSE and the residuals RES.
%
%         [R,XYZ1,XYZ2,SSE,RES,EXITFLAG] = CYL_FIT(RI,XYZ1I,XYZ2I,XYZP)
%         returns the EXITFLAG from MATLAB function lsqnonlin or
%         fminsearch.  Use "help lsqnonlin" or "help fminsearch" to see
%         the exit conditions for the different values of EXITFLAG.
%
%         NOTES:  1.  See cyl_plt.m to plot the resulting cylinder.
%
%                 2.  The M-file dist2cyl.m and pts2lin.m must be in
%                 the current path or directory.
%
%         17-Jul-2013 * Mack Gardner-Morse
%
%         21-Nov-2014 * Mack Gardner-Morse * Added the METHOD option.
%
%         25-Nov-2014 * Mack Gardner-Morse * Check cylinder axis matches
%                                            input axis direction.
%
%         10-Aug-2022 * Mack Gardner-Morse * Removed the METHOD option.
%         Only uses the Matlab nonlinear least squares function
%         lsqnonlin.
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  error(' *** ERROR in CYL_FIT:  Four (4) inputs are required!');
end
%
if (nargin<5)||isempty(tol)
  tol = eps^(2/3);
end
%
% if (nargin<6)||isempty(method)
%   method = true;
% end
%
% Check Inputs
%
[nrow1,ncol1] = size(xyz1i);
[nrow2,ncol2] = size(xyz2i);
[~,ncol3] = size(xyzp);
if nrow1~=1||nrow2~=1||ncol1~=3||ncol2~=3||ncol3~=3
  error([' *** ERROR in CYL_FIT:  All input points must have', ...
         ' three columns!']);
end
%
% Get Cylinder Parameters
%
xi = zeros(7,1);
if ri>0
  xi(1) = ri;
else
  xi(1) = 1;
end
xi(2:4) = xyz1i';
xi(5:7) = xyz2i';
zci = xyz2i-xyz1i;                     % Direction of the cylinder axis
zci = zci./sqrt(zci*zci');             % Unit direction vector
%
% Set Optimization Parameters
%
% if method
%   opt = optimset('fminsearch');
% else
opt = optimset('lsqnonlin');
opt = optimset(opt,'Display','off','Algorithm','levenberg-marquardt');
%   opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
%                  'MaxIter',2e+3,'Algorithm','levenberg-marquardt');
% end
%
opt.MaxFunEvals = 5e+7;
opt.MaxIter = 2e+7;
opt.TolFun = tol;
opt.TolX = tol;
%
% Fit Cylinder
%
% if method
%   [x,~,exitflag] = fminsearch(@(x) sum(dist2cyl(x,xyzp).^2),xi,opt);
%   res = dist2cyl(x,xyzp);
%   sse = res'*res;
% else
lb = -Inf(7,1);
lb(1) = 0;              % Radius must be positive
[x,sse,res,exitflag] = lsqnonlin(@(x) dist2cyl(x,xyzp),xi,lb,[],opt);
% end
%
r = x(1);
xyz1 = x(2:4)';
xyz2 = x(5:7)';
%
% Check Cylinder Direction Matches Input Direction
%
zc = xyz2-xyz1;                        % Direction of the cylinder axis
zc = zc./sqrt(zc*zc');                 % Unit direction vector
cdir = zci*zc';                        % Positive == match input direction
if cdir<0
  xyzt = xyz1;
  xyz1 = xyz2;
  xyz2 = xyzt;
  clear xyzt;
  zc = -zc;
end
%
% Get Ends of the Cylinder that Fit the Range of the Data
%
[xyzl,t] = pts2lin(xyz1,zc,xyzp);      % Get points and distances along axis
[~,idmn] = min(t);                     % t for first end of cylinder
[~,idmx] = max(t);                     % t for second end of cylinder
xyz1 = xyzl(idmn,:);                   % Coordinates of first end of cylinder
xyz2 = xyzl(idmx,:);                   % Coordinates of second end of cylinder
%
return