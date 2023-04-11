function dats = comb_dat(datl,datm,datt);
%COMB_DAT Combines the sagittal lateral and medial condyle and trochlea
%         MRI knee data cell arrays into a single data cell array.
%
%         DATS = COMB_DAT(DATL,DATM,DATT) given the data cell arrays
%         containing three (3) columns matrices with slice coordinate
%         point data for the lateral condyle, DATL, medial condyle,
%         DATM, and trochlea, DATT, returns the cell array, DATS, that
%         contains the combined ordered slice data from the three input
%         data cell arrays.
%
%         NOTES:  1.  To insure the correct ordering of the slices,
%                 the input cell arrays must be input in the correct
%                 order (lateral, medial and trochlea).
%
%                 2.  The function assumes all the slices are in order
%                 within the data cell arrays.
%
%                 3.  The M-files plane_fit.m and sl_dir.m must be in
%                 the current path or directory.
%
%         03-Jul-2014 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in COMB_DAT:  No or not enough input data!');
end
%
% Get Direction Between the Centers of the Condyles
%
xyzl = cell2mat(datl);
xyzm = cell2mat(datm);
%
vx = mean(xyzl)-mean(xyzm);
%
% Get Directions Within the Data Cell Arrays
%
slvm = sl_dir(datm);
slvt = sl_dir(datt);
slvl = sl_dir(datl);
%
% Reverse Slice Order if Necessary
%
if slvm*vx'<0
  datm = flipud(datm);
end
%
if slvt*vx'<0
  datt = flipud(datt);
end
%
if slvl*vx'<0
  datl = flipud(datl);
end
%
% Combine Slice Data (Cell Arrays)
%
dats = [datm; datt; datl];
%
return