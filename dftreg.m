function [reg_img,drow,dcol,err,phase] = dftreg(fix_img,mov_img,prec)
%DFTREG    Image registration by crosscorrelation for in-plane 
%          translations.
%
%          REG_IMG = DFTREG(FIX_IMG,MOV_IMG) given the 2D reference
%          image matrix, FIX_IMG, and the image to be transformed,
%          MOV_IMG, the registered (transformed) image, REG_IMG is
%          returned.
%
%          [REG_IMG,DROW,DCOL,ERR,PHASE] = DFTREG(FIX_IMG,MOV_IMG)
%          returns the row (image X, plot Y) translation, DROW, the
%          the column (image Y, plot X) translation, DCOL, the
%          translation invariant normalized RMS error between the
%          reference and registered images, ERR, and the phase between
%          the reference and moving images, PHASE (usually zero).
%
%          NOTES:  1.  This is a wrapper function for dftregistration.m
%                  which is a Mathworks File Exchange function.  This
%                  function does the 2D FFTs and inverse 2D FFT on the
%                  input images and output image.
%
%                  2.  If the input images are not double precision,
%                  the images are converted to double using the
%                  Matlab function im2double. 
%
%                  3.  The Matlab file dftregistration.m must be in the
%                  current directory or path.
%
%          30-Jun-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error([' *** ERROR in dftreg:  Two input image matrices are', ...
         ' required!']);
end
%
if (nargin<3)
  prec = 100;           % Default precision of translations (1/100)
end
%
% Make Input Image Matrices Double Precision
%
fix_img = im2double(fix_img);
mov_img = im2double(mov_img);
%
% Call dftregistration.m
%
[out,reg_img] = dftregistration(fft2(fix_img),fft2(mov_img),prec);
%
% Parse Output
%
reg_img = abs(ifft2(reg_img));
%
drow = out(3);
dcol = out(4);
err = out(1);
phase = out(2);
%
return