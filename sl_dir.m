function [slvec,slx] = sl_dir(dat,iall);
%SL_DIR   Gets the slice direction between MRI slices in the ordered 
%         slice data from the digitized MRI knee cell arrays.
%
%         SLVEC = SL_DIR(DAT) given a cell array containing three (3)
%         columns matrices with slice coordinate point data, DAT,
%         returns the three (3) column slice separation vector, SLVEC,
%         between the second and first slice.
%
%         [SLVEC,SLX] = SL_DIR(DAT) returns the magnitude (separation
%         distance) of the separation vector, SLX.
%
%         [SLVEC,SLX] = SL_DIR(DAT,IALL) if the logical scalar, IALL, is
%         true, returns the separation vector, SLVEC, between each slice
%         and the separation distance, SLX between each slice in matrix
%         with the different slice data in rows.
%
%         NOTES:  1.  The M-file plane_fit.m must be in the current path
%                 or directory.
%
%                 2.  Only checks the direction between the first two
%                 slices.
%
%         03-Jul-2014 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in SL_DIR:  No input data!');
end
%
if (nargin<2)
  iall = false;
end
%
dat = dat(:);
if iall
  nslice = size(dat,1);
else
  nslice = 2;
end
cntr = zeros(nslice,3);
nvec = zeros(3,nslice);
%
% Fit Planes to Slice Data
%
for k = 1:nslice
   xyz = dat{k};
   [cntr(k,:),nvec(:,k)] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));
   vec = xyz(end,:)-xyz(1,:);
   vec1(k,:) = vec./norm(vec);
end
%
slx = zeros(nslice-1,1);
slvec = zeros(nslice-1,3);
%
for k = 2:nslice
%
% Slice Separations
%
   slx(k-1) = (cntr(k,:)-cntr(k-1,:))*nvec(:,k-1);
   slvec(k-1,:) = slx(k-1)*nvec(:,k-1)';
%
end
%
slx = abs(slx);         % Scalar distance (always positive)
%
return