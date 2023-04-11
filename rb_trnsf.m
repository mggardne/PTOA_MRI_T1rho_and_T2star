function [fr,fr3,dr] = rb_trnsf(s,r,t,f3fc1,f3fb1,f3fc2,f3fb2)
%RB_TRNSF  Given the parameters for a rigid body transformation with
%          scaling, calculates the transformation for two pairs of
%          knee cartilage two-dimensional (2D) or three-dimensional (3D)
%          coordinates which are grouped and returned in a cell array.
%
%          FR = RB_TRNSF(S,R,T,F3FC1,F3FB1,F3FC2,F3FB2) Inputs are the
%          scaling, S, rotation matrix , R, translation vector, T, and
%          two pairs of knee cartilage coordinates.  The first pair is
%          the cartilage surface, F3FC1, and cartilage-bone interface,
%          F3FB1, coordinates.  The second pair is the cartilage
%          surface, F3FC2, and cartilage-bone interface, F3FB2,
%          coordinates, returns a cell array, FR, with the transformed
%          coordinates.  The cartilage coordinates are in the first
%          row of FR, the first pair of coordinates are in the first
%          column, and similarly, the second pair of coordinates are in
%          the second column.
%
%          [FR,FR3,DR] = RB_TRNSF(S,R,T,F3FC1,F3FB1,F3FC2,F3FB2)
%          Returns the 3D coordinates in a cell array, FR3 (similar to
%          FR), and a three-dimensional matrix with the unit vectors of
%          the directions of the coordinates (last point minus the
%          first point), DR.  The first row is cartilage and the second
%          row is bone.  The columns are the three vector components,
%          and the third index is for the pairs of coordinates.
%
%          NOTES:  1.  As long as the transformation parameters' sizes
%                  match the dimensions of the coordinates, the
%                  function works for 2D or 3D coordinates.
%
%          13-Jan-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<7)
  error([' *** ERROR in RB_TRNSF:  Seven (7) input variables are', ...
         ' required!']);
end
%
% Initialize the Cell Arrays and Direction Vector Matrix
%
fr3 = cell(2,2);
dr = zeros(2,3,2);      % Direction vectors - 1st column (1 - cartilage, 2 - bone)
fr = cell(2,2);
%
fr3{2,1} = f3fb1;
nf3d = size(f3fb1,1);   % Number of points
dr(2,:,1) = f3fb1(nf3d,:)-f3fb1(1,:);
fr{2,1} = s*f3fb1*r+repmat(t,nf3d,1);  % Transform to pixel coordinates
fr{2,1} = fr{2,1}(:,1:2);              % Convert from 3D to 2D
%
fr3{2,2} = f3fb2;
nf3d = size(f3fb2,1);   % Number of points
dr(2,:,2) = f3fb2(nf3d,:)-f3fb2(1,:);
fr{2,2} = s*f3fb2*r+repmat(t,nf3d,1);  % Transform to pixel coordinates
fr{2,2} = fr{2,2}(:,1:2);              % Convert from 3D to 2D
%
fr3{1,1} = f3fc1;
nf3d = size(f3fc1,1);   % Number of points
dr(1,:,1) = f3fc1(nf3d,:)-f3fc1(1,:);
fr{1,1} = s*f3fc1*r+repmat(t,nf3d,1);  % Transform to pixel coordinates
fr{1,1} = fr{1,1}(:,1:2);              % Convert from 3D to 2D
%
fr3{1,2} = f3fc2;
nf3d = size(f3fc2,1);   % Number of points
dr(1,:,2) = f3fc2(nf3d,:)-f3fc2(1,:);
fr{1,2} = s*f3fc2*r+repmat(t,nf3d,1);  % Transform to pixel coordinates
fr{1,2} = fr{1,2}(:,1:2);              % Convert from 3D to 2D
%
return