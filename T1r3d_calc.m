function [T1r,T1r_amp,sse,exit_flags] = T1r3d_calc(T1I,slt,T1rho0)
%T1R3D_CALC Given a three-dimensional matrix of T1/T2 intensities from
%          a MRI image slice and a vector of the spin lock/echo times, 
%          calculates the T1rho/T2* map for the slice.
%
%          T1R = T1R3D_CALC(T1I,SLT) Given a three-dimensional matrix
%          of T1/T2 intensities from a MRI image slice, T1I, and a
%          vector of the spin lock/echo times, SLT, calculates the
%          T1rho/T2* map for the slice, T1R.  The matrix of T1/T2
%          intensities, T1I, must be a three dimensional array of the
%          image slices acquired at different spin lock/echo times.
%          The third dimension being the index to the intensities from
%          different spin lock/echo times.
%
%          T1R = T1R3D_CALC(T1I,SLT,T1rho0) Uses the value of T1rho0 as
%          the initial T1rho/T2* value for fitting the monoexponential
%          if the weighted linear least squares values are less than
%          0.1 or greater than T1rho0.  The default initial T1rho/T2*
%          value, T1rho0, is 80.  This is near the maximum value for
%          T1rho in human knee cartilage.
%
%          [T1R,T1R_AMP,SSE,EXIT_FLAGS] = T1R3D_CALC(T1I,SLT) Returns
%          a matrix of the amplitude of the exponential fits, T1R_AMP,
%          a matrix of sum of squared errors of the fits, SSE, and
%          a matrix of the Matlab exit flag from the nonlinear least
%          squares algorithm, EXIT_FLAGS.
%
%          NOTES:  1.  The T1rho/T2* map is calculated as a mono-
%                  exponential as a function of spin lock/echo times.
%
%                  2.  The images slice data and spin lock/echo times
%                  are assumed to be in ascending order with increasing
%                  index into the third dimension.
%
%                  3.  T1/T2 pixels with values less than 0.1 in the
%                  image from the first spin lock/echo time are not fit
%                  and set to zero.
%
%                  4.  M-file exp_fun1.m must be in the current
%                  directory or path.
%
%                  5.  The Matlab Optimization toolbox is required.
%                  The nonlinear least squares is performed by the
%                  Matlab function lsqcurvefit in the optimization
%                  toolbox.
%
%                  6.  The nonlinear least squares uses the Levenberg-
%                  Marquardt algorithm in the Matlab function
%                  lsqcurvefit.
%
%                  7.  The Matlab Parallel Computing toolbox is 
%                  required.  The Matlab parallel construct parfor is
%                  used to calculate the T1rho/T2* values in parallel.
%                  Use Matlab command parpool to control the number of 
%                  workers.
%
%          28-Apr-2021 * Mack Gardner-Morse
%
%          05-Jul-2021 * Mack Gardner-Morse * Generalized comments,
%          eliminated "nd" variable and corrected linear least squares
%          amplitude estimate by taking exponential (exp function).
%

%#######################################################################
%
% Check for Inputs
%
if nargin<2
  error([' *** ERROR in T1r3d_calc:  A vector of echo times is', ...
         ' required!']);
end
%
if nargin<3
  T1rho0 = 80;          % Initial T1rho/T2* value (maximum expected)
end
%
if isempty(T1rho0)
  T1rho0 = 80;          % Initial T1rho/T2* value (maximum expected)
end
%
nslt = size(T1I,3);     % Number of spin lock/echo times
%
if nslt==1
  error([' *** ERROR in T1r3d_calc:  Input T1/T2 intensity matrix', ...
         ' must have three dimensions!']);
end
%
% Get Image Information
%
isz = size(T1I,1:2);    % Size of image
npx = prod(isz);        % Number of pixels in the image 
%
slt = double(slt(:));   % Spin lock times as a column vector
%
if size(slt,1)~=nslt
  error([' *** ERROR in T1r3d_calc:  The size of the third', ...
         ' dimension of the T1/T2 intensity matrix must match the', ...
         ' number of echo times!']);
end
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter', ...
               2e+3,'Algorithm','levenberg-marquardt','Jacobian', ...
               'on','UseParallel',true);
%
fun = @exp_fun1;        % Exponential function
%
% Don't Fit Pixels with Small (Zero) Initial Values
%
dat3d = double(T1I);    % Do calculations in double precision
vmat = dat3d(:,:,1)>0.1;     % > small value for initial spin lock time
vmat = repmat(vmat,1,1,nslt);          % Valid values of intensities
vmat = reshape(vmat,npx,nslt)';
%
% Set Up Arrays for Loop
%
xdat = [slt ones(nslt,1)];             % Spin lock/echo times and exponential amplitude
dat3d = reshape(dat3d,npx,nslt)';      % Intensities as a vector
%
T1r = zeros(npx,1,'single');           % Nonlinear least squares time constants
T1r_amp = zeros(npx,1,'single');       % Nonlinear least squares amplitudes
sse = zeros(npx,1,'single');           % Nonlinear least squares sum of squared errors
exit_flags = zeros(npx,1,'int8');      % Nonlinear least squares exit flags
%
% Get Log for Linear Least Squares
%
dat3dl = double(vmat);
dat3dl(vmat) = log(dat3d(vmat));
%
% Calculate Linear Least Squares T1rho/T2* Values
%
p = xdat\dat3dl;
p = abs(p);             % Get magnitude of complex numbers
T1r_lin = -1./p(1,:);   % Time constant in ms
lin_amp = exp(p(2,:));  % Amplitude
%
% Nonlinear Least Squares Exponential Fit to Get T1rho/T2* Values
%
parfor k = 1:npx        % Pixels
   if vmat(1,k)         % Only calculate if valid data
     ydat = dat3d(:,k);
     if T1r_lin(k)>0.1&&T1r_lin(k)<=T1rho0  % Use linear least squares as initial parameters
       rp0 = [lin_amp(k) T1r_lin(k)];
     else
       rp0 = [max(ydat) T1rho0];
     end
%
     [rp,~,~,eflag] = lsqcurvefit(fun,rp0,slt,ydat,[],[],opt);
%
     if isinf(rp(2))
       T1r(k) = single(0);   % Set Inf to zero
     else
       T1r(k) = single(rp(2));
     end
     T1r_amp(k) = single(rp(1));
     d = exp_fun1(rp,slt)-ydat;
     sse(k) = single(sqrt(d'*d));
     exit_flags(k) = int8(eflag);
   end
end
%
% Reshape Vectors into Matrices
%
T1r = reshape(T1r,isz);
T1r_amp = reshape(T1r_amp,isz);
sse = reshape(sse,isz);
exit_flags = reshape(exit_flags,isz);
%
return