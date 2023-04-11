function [tc,amp,rss,npx,id,tcpx,amppx,rsspx,nps] = psl_ana(v, ...
              masklay,maskroi,rsl,nrsl,rslbs,time,ntime,fun,init,t0,opt)
%PSL_ANA   Given a four-dimensional matrix of T1/T2 intensities from a
%          MRI image volume, a mask, and a vector of spin lock/echo
%          times, calculates the T1rho/T2* map for the volume.
%
%          TC = PSL_ANA(V,MASKLAY,MASKROI,RSL,NRSL,RSLBS,TIME,NTIME,
%          FUN,INIT) Given a four-dimensional matrix of T1/T2
%          intensities from a MRI image volume, V, where the first two
%          dimensions are an image, the third dimension are the slices
%          and the fourth dimension are the spin lock/echo times, three-
%          dimensional logical masks with the first dimension being the
%          image, the second dimension being the superficial layer in
%          the first column and deep layer in the second column and the
%          third dimension being slices in a cell array of masks with
%          the cell index to the femur and tibia, MASKLAY, two-
%          dimensional logical masks with the first dimension being the
%          image, the second dimension being slices in a cell array of
%          masks with the first index to the femur and tibia and the
%          second index to lateral and medial, MASKROI, array of slices
%          with segmentations, RSL, the number of slices with
%          segmentations, NRSLS, cell array with slices in the femur
%          and tibia, RSLBS, vector of spin lock/echo times, TIME, the
%          number of spin lock/echo times, NTIME, a function handle,
%          FUN, to evaluate the monoexponential and its derivatives,
%          and initialization flag, INIT, to determine the initial
%          starting parameters, calculates T1rho/T2* for the two
%          compartments (lateral and medial), two bones (femur and
%          tibia), three regions of interest (ROIs) (anterior/trochlea,
%          central, and posterior), and two layers (superficial and
%          deep) in array TC.
%
%          TC = PSL_ANA(V,MASKLAY,MASKROI,RSL,NRSL,RSLBS,TIME,NTIME,
%          FUN,INIT,T0) If INIT is greater than zero, uses T0 as the
%          initial default time constant for the nonlinear exponential
%          fit.  The default value is 50.  The mean value of the
%          maximum intensities are used as the initial amplitude.  If
%          INIT is equal to zero, the exponential equation is
%          linearized using logarithms and least squares are used to
%          solve for the initial parameters for the nonlinear least
%          squares.  If INIT is less than zero, weighted least squares
%          using the intensities is used to solve for the starting
%          parameters.
%
%          TC = PSL_ANA(V,MASKLAY,MASKROI,RSL,NRSL,RSLBS,TIME,NTIME,
%          NTIME,FUN,INIT,T0,OPT) a structure, OPT, with options for
%          the curvefit solver.  See Matlab command "optimset" for a
%          list of parameters.
%
%          [TC,AMP,RSS,NPX,ID] = PSL_ANA(...) returns the amplitudes
%          of the exponential fits, AMP, the sum of squared residuals,
%          RSS, the number of pixels in the curvefits, NPX, and a four-
%          dimensional array with the first dimension being the
%          compartment (0 - lateral and 1 - medial), the second
%          dimension being bone (0 - femur and 1 - tibia), the third
%          dimension being ROIs (0 - anterior/trochlea, 1 - central,
%          and 2 - posterior), and the fourth dimension being layer
%          (0 - deep and 1 - superficial) that identifies the analyzes,
%          ID. 
%
%          [TC,AMP,RSS,NPX,ID,TCPX,AMPPX,RSSPX,NPS] = PSL_ANA(...)
%          returns the T1rho/T2* for each pixel in the two compartments
%          (lateral and medial), two bones (femur and tibia), three ROIs
%          (anterior/trochlea, central, and posterior) and two layers
%          (superficial and deep) in a cell array, TCPX, the
%          exponential amplitudes for each pixel in a cell array, AMPPX,
%          the sum of squared residuals for each pixel, RSSPX, and the
%          number of pixels in each slice within the regions in a cell
%          array, NPS.
%
%          NOTES:  1.  The T1rho/T2* map is calculated as a mono-
%                  exponential as a function of spin lock/echo times.
%
%                  2.  The images volume data and spin lock/echo times
%                  are assumed to be in ascending order with increasing
%                  index into the fourth dimension of the image data V.
%
%                  3.  The function handle, FUN, is to the M-file
%                  exp_fun1.m which must be in the current directory or
%                  path.
%
%                  4.  T0 must be between 1 and 100.
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
%          20-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<10
  error([' *** ERROR in psl_ana:  Ten input variables are', ...
         ' required!']);
end
%
if nargin<11
  t0 = 50;              % Initial estimate of time constant
end
%
% Check Initial Time Constant
%
if isempty(t0)||t0<1||t0>100
  t0 = 50;
end
%
% Default Optimization Parameters
%
if nargin<12
  opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
                 'MaxIter',2e+3,'Algorithm','levenberg-marquardt', ...
                 'Jacobian','on','UseParallel',true);
end
%
% Initialize Variables
%
% Columns in variable "id":
%   1 - Bones:  0 - Femur and 1 - Tibia
%   2 - Compartments:  0 - Lateral and 1 - Medial
%   3 - ROIs:  0 - Anterior/Trochlea, 1 - Central and 2 - Posterior
%   4 - Layers:  0 - Deep and 1 - Superficial
%
tc = zeros(24,1);       % Time constant - T1rho/T2*
amp = zeros(24,1);      % Amplitude of exponential
rss = zeros(24,1);      % Sum of squared residuals
npx = zeros(24,1);      % Number of pixels in regional fit
id = zeros(24,4);       % Columns: 1 - Bone, 2 - Compartment, 3 - ROI and 4 - Layer
%
nps = cell(24,1);       % Number of pixels in each slice
tcpx = cell(24,1);      % Time constant - T1rho/T2* for each pixel
amppx = cell(24,1);     % Amplitude of exponential for each pixel
rsspx = cell(24,1);     % Sum of squared residuals for each pixel
%
% Loop through Slices
% rimgs dimensions: 1 - Slice, 2 - Bone, 3 - Compartment, 4 - ROI, and
%                   5 - Layer
%
rimgs = cell(nrsl,2,2,3,2);
np = zeros(2,2,3,2,nrsl);
%
for k = 1:nrsl
%
   ksl = rsl(k);        % Slice number
%
   rimg = squeeze(v(:,:,ksl,:));  % T1/T2 data for slice and spin lock/echo time
%
% Loop through Bones
%
   for kb = 1:2         % 1 - Femur and 2 - Tibia
%
      mskb = maskroi{kb};    % ROI masks for this bone
      rslb = rslbs{kb};      % Slices for this bone
      idxb = rslb==ksl;      % Index to slice in rslb/mskb
      if ~any(idxb)
        continue;       % No bone segmentation on this slice
      end
%
      maskb = masklay{kb};        % Layer mask for this bone
%
% Loop through Compartments
%
      for kc = 1:2      % 1 - Lateral and 2 - Medial
%
         mskc = mskb{kc};    % ROI masks for this compartment
%
% Loop through ROIs
%
         for kr = 1:3   % 1 - anterior/trochlea, 2 - central and 3 - posterior
%
            mskr = mskc{kr};
            msks = mskr(:,idxb);
%
            npps = sum(msks);
            if npps==0
              continue; % No ROI on this slice
            end
%
% Loop through Layers
%
            for kl = 1:2     % 1 - superficial, 2 - deep
%
               maskl = maskb(:,kl,k);
               mskl = msks&maskl;
%
               npps = sum(mskl);
               np(kb,kc,kr,kl,k) = npps;
%
               rimgl = rimg(repmat(mskl,1,1,ntime));
               rimgs{k,kb,kc,kr,kl} = reshape(rimgl,npps,ntime)';
%
            end         % End of kl loop - layers loop
%
         end            % End of kr loop - ROIs loop
%
      end               % End of kc loop - compartments loop
%
   end                  % End of kb loop - bones loop
%
end                     % End of k loop - slices loop
%
clear npps rimg rimgl;
%
% Calculate T1rho/T2* by Looping through Bones, Compartments, ROIs, and
% Layers
%
for kb = 1:2            % 1 - Femur and 2 - Tibia
%
   ikb = kb-1;          % Bone code identifier (0 - femur, 1 - tibia)
%
   for kc = 1:2         % 1 - Lateral and 2 - Medial
%
      ikc = kc-1;       % Compartment code (0 - lateral, 1 - medial)
%
      for kr = 1:3      % 1 - anterior/trochlea, 2 - central and 3 - posterior
%
         ikr = kr-1;    % 0 - anterior/trochlea, 1 - central and 2 - posterior
%
         for kl = 1:2   % 1 - superficial, 2 - deep
%
            ikl = 2-kl; % 0 - deep, 1 -superficial
%
            idx = 12*kb+6*kc+2*kr+kl-20;    % Index to variables
            nps{idx} = squeeze(np(kb,kc,kr,kl,:));    % Number of pixels in slice
%
            rimgl = cat(2,rimgs{:,kb,kc,kr,kl});
            npx(idx) = size(rimgl,2);
            r = rimgl(:);
%
            if npx(idx)==0
              continue; % No pixels!
            end
%
            amp0 = mean(rimgl(1,:));
            rp0c = [amp0; t0];         % Constant initial parameters
%
            tl = repmat(time,npx(idx),1);
%
            if init==0  % Use linear least squares for initial constants
              xdat = [tl ones(ntime*npx(idx),1)];
              lr = log(r);
              rpl = xdat\lr;
              rpl = abs(rpl);     % Get magnitude of complex numbers
              t0l = -1/rpl(1);
              if t0l>1&&t0l<100
                rp0 = [exp(rpl(2)); t0l];
              else
                rp0 = rp0c;
              end
            elseif init<0    % Use weighted least squares for initial constants
              xdat = [tl.*r r];
              lr = log(r);
              rpl = xdat\(lr.*r);
              rpl = abs(rpl);     % Get magnitude of complex numbers
              t0l = -1/rpl(1);
              if t0l>1&&t0l<100
                rp0 = [exp(rpl(2)); t0l];
              else
                rp0 = rp0c;
              end
            else        % Use constant initial parameters
              rp0 = rp0c;
            end
%
% Nonlinear Least Squares Exponential Fit to Get T1rho/T2* Value
%
            [rp,err] = lsqcurvefit(fun,rp0,tl,r,[],[],opt);
%
            tc(idx) = rp(2);           % Time constant - T1rho/T2* 
            amp(idx) = rp(1);          % Amplitude of exponential
            rss(idx) = err;            % Residual sum of squares
            id(idx,:) = [ikb ikc ikr ikl];  % Compartment, bone, ROI, and layer IDs
%
% Calculate T1rho/T2* for Each Pixel in Each Layer
%
            if nargout>5
              tcp = zeros(npx(idx),1);
              ampp = tcp;
              rssp = tcp;
%
% Get Initial Parameters
%
              if init==0     % Use linear least squares for initial constants
                xdat = [time ones(ntime,1)];
                rpa = xdat\log(rimgl);
                rpa = abs(rpa);        % Get magnitude of complex numbers
                t0a = -1./rpa(1,:);
                ampa = exp(rpa(2,:));
                rpa = [ampa; t0a];
                idv = t0a<1|t0a>100;
                if any(idv)
                  rpa(:,idv) = repmat(rp0c,1,sum(idv));
                end
              elseif init<0  % Use weighted least squares for initial constants
                rpa = repmat(rp0c,1,npx(idx));
                for p = 1:npx(idx)
                   rimgp = rimgl(:,p);
                   xdat = [time.*rimgp rimgp];
                   lr = log(rimgp);
                   rpl = xdat\(lr.*rimgp);
                   rpl = abs(rpl);     % Get magnitude of complex numbers
                   t0l = -1/rpl(1);
                   if t0l>=1&&t0l<=100
                     rpa(:,p) = [exp(rpl(2)); t0l];
                   end
                end
              else           % Use constant initial parameters
                rpa = repmat(rp0c,1,npx(idx));
              end
%
% Nonlinear Least Squares Exponential Fit to Get T1rho/T2* Pixel Values
%
              parfor p = 1:npx(idx)
%      
                 rimgp = rimgl(:,p);
%      
                 [rp,err] = lsqcurvefit(fun,rpa(:,p),time,rimgp,[],[],opt);
                 tcp(p) = rp(2);
                 ampp(p) = rp(1);
                 rssp(p) = err;
%
              end
%
% Save Results to Outputs
%
              tcpx{idx} = tcp;
              amppx{idx} = ampp;
              rsspx{idx} = rssp;
%
            end
%
         end            % End of kl loop - layers loop
%
      end               % End of kr loop - ROIs loop
%
   end                  % End of kc loop - compartments loop
%
end                     % End of kb loop - bones loop
%
return