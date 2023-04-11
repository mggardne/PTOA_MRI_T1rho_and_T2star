function [varargout] = coord_tf(xyz0,rotmat,varargin)
%COORD_TF Transforms MRI knee data cell arrays by a coordinate 
%         transformation (translation and rotation).
%
%         DATT = COORD_TF(XYZ0,ROTMAT,DAT) given the axes offset, XYZ0,
%         and rotation matrix, ROTMAT, transforms the coordinates in
%         MRI knee data cell array, DAT, by translating and rotating
%         the coordinates (xyzt = (xyz-xyz0)*rotmat).  The transformed
%         coordinates are returned in a cell array, DATT.
%
%         [DAT1T,DAT2T,DAT3T,...] = COORD_TF(XYZ0,ROTMAT,DAT1,DAT2,DAT3,
%         ...) returns additional cell arrays that have been transformed
%         by the same coordinate transformation.
%
%         NOTES:  1.  Inputs must be a MRI knee data cell arrays with
%                 the coordinate matrices for each slice in rows in the
%                 cell array.
%
%                 2.  See tf_coord.m for a transformation that first
%                 rotates and then translates the coordinates.
%
%         24-Jul-2014 * Mack Gardner-Morse
%
%         28-Sep-2022 * Mack Gardner-Morse * Replaced loop and function
%                       call with cell array functions cellfun and
%                       mat2cell.
%

%#######################################################################
%
% Check Inputs
%
if nargin<3
  error([' *** ERROR in COORD_TF:  At least one input MRI knee', ...
         ' data cell array is required!']);
end
%
xyz0 = xyz0(:)';
if size(xyz0,2)~=3
  error(' *** ERROR in COORD_TF:  XYZ0 must be of length three (3)!');
end
%
if size(rotmat,1)~=3||size(rotmat,2)~=3
  error([' *** ERROR in COORD_TF:  ROTMAT must be a 3x3 rotation', ...
         '  matrix!']);
end
%
iloop = nargin-2;
for k = 1:iloop
   dat = varargin{k};
   if ~iscell(dat)
     error(' *** ERROR in COORD_TF:  Inputs must be cell arrays!');
   end
end
%
% Check Outputs
%
if nargout~=nargin-2
  error([' *** ERROR in COORD_TF:  Number of inputs do not match', ...
        ' the number of outputs!']);
end
%
% Transform Cell Arrays
%
for k = iloop:-1:1
%
% Get Cell Array
%
   dat = varargin{k};
   nsiz = size(dat);
   dat = dat(:);
%
% Get Slice Information
%
   nsd = cellfun('size',dat,1);
%
% Get and Rotate All of the Coordinates
%
   xyz = cell2mat(dat);
%
   nd = size(xyz,1);
%
   toff = repmat(xyz0,nd,1);
   txyz = (xyz-toff)*rotmat;
%
% Reform Cell Arrays
%
   datt = mat2cell(txyz,nsd);
   datt = reshape(datt,nsiz);
%
% Output Transformed Coordinates
%   
   varargout{k} = datt;
%
end
%
return