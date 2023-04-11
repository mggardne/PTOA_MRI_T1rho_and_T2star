function [s,r,t] = decomp(x,y,iscale)
%
% See Challis JH: A procedure for determining rigid body
% transformation parameters.  J Biomech 28(6):733-7,1995.
%
if nargin<3
  iscale = true;
end
%
[m,n] = size(x);
o = ones(m,1);
mx = mean(x);
my = mean(y);
xc = x-o*mx;
yc = y-o*my;
%
[p,d,q] = svd(xc'*yc);
%
r = p*q';
if (det(r)<0)
   p(:,3) = -p(:,3);
   r = p*q';
end
%
if iscale
  s = trace(d)/trace(xc'*xc);
else
  s = 1;
end
%
t = my'-s*r'*mx';
t = t';
%
return