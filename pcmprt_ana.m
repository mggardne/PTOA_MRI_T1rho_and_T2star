function [tc,amp,rss,npx,id,tcpx,amppx,rsspx,nps] = pcmprt_ana(v, ...
        masklay,maskroi,rsls,nrsls,rsl,rslbs,time,ntime,fun,init,t0,opt)
%PCMPRT_ANA Given a four-dimensional matrix of T1/T2 intensities from
%          a MRI image volume, a mask, and a vector of spin lock/echo
%          times, calculates the T1rho/T2* map for the volume.
%
%          TC = PCMPRT_ANA(V,MASKLAY,MASKROI,RSLS,NRSLS,RSL,RSLBS,TIME,
%          NTIME,FUN,INIT) Given a four-dimensional matrix of T1/T2
%          intensities from a MRI image volume, V, where the first two
%          dimensions are an image, the third dimension are the slices
%          and the fourth dimension are the spin lock/echo times, three-
%          dimensional logical masks with the first dimension being the
%          image, the second dimension being the superficial layer in
%          the first column and deep layer in the second column and the
%          third dimension being slices in a cell array of masks with
%          the cell index to femur and tibia, MASKLAY, two-dimensional
%          logical masks with the first dimension being the image, the
%          second dimension being slices in a cell array of masks with
%          the first index to lateral and medial and the second index
%          to the femur and tibia, MASKROI, cell array with the slices
%          within each compartment, RSLS, the number of slices in each
%          compartment, NRSLS, array of slices with segmentations, RSL,
%          cell array with slices in the femur and tibia, RSLBS, vector
%          of spin lock/echo times, TIME, the number of spin lock/echo
%          times, NTIME, a function handle, FUN, to evaluate the
%          monoexponential and its derivatives, and initialization flag,
%          INIT, to determine the initial starting parameters,
%          calculates T1rho/T2* for the two compartments (lateral and
%          medial), two bones (femur and tibia), three regions of
%          interest (ROIs) (anterior/trochlea, central, and posterior),
%          and two layers (superficial and deep) in array TC.
%
%          TC = PCMPRT_ANA(V,MASKLAY,MASKROI,RSLS,NRSLS,RSL,RSLBS,TIME,
%          NTIME,FUN,INIT,T0) If INIT is greater than zero, uses T0 as
%          the initial default time constant for the nonlinear
%          exponential fit.  The default value is 50.  The mean value
%          of the maximum intensities are used as the initial amplitude.
%          If INIT is equal to zero, the exponential equation is
%          linearized using logarithms and least squares are used to
%          solve for the initial parameters for the nonlinear least
%          squares.  If INIT is less than zero, weighted least squares
%          using the intensities is used to solve for the starting
%          parameters.
%
%          TC = PCMPRT_ANA(V,MASKLAY,MASKROI,RSLS,NRSLS,RSL,RSLBS,TIME,
%          NTIME,FUN,INIT,T0,OPT) a structure, OPT, with options for
%          the curvefit solver.  See Matlab command "optimset" for a
%          list of parameters.
%
%          [TC,AMP,RSS,NPX,ID] = PCMPRT_ANA(...) returns the amplitudes
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
%          [TC,AMP,RSS,NPX,ID,TCPX,AMPPX,RSSPX,NPS] = PCMPRT_ANA(...)
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
%          10-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<11
  error([' *** ERROR in pcmprt_ana:  Eleven input variables are', ...
         ' required!']);
end
%
if nargin<12
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
if nargin<13
  opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
                 'MaxIter',2e+3,'Algorithm','levenberg-marquardt', ...
                 'Jacobian','on','UseParallel',true);
end
%
% Initialize Variables
%
% Columns in variable "id":
%   1 - Compartments:  0 - Lateral and 1 - Medial
%   2 - Bones:  0 - Femur and 1 - Tibia
%   3 - ROIs:  0 - Anterior/Trochlea, 1 - Central and 2 - Posterior
%   4 - Layers:  0 - Deep and 1 - Superficial
%
tc = zeros(24,1);       % Time constant - T1rho/T2*
amp = zeros(24,1);      % Amplitude of exponential
rss = zeros(24,1);      % Sum of squared residuals
npx = zeros(24,1);      % Number of pixels in regional fit
id = zeros(24,4);       % Columns: 1 - Compartment, 2 - Bone, 3 - ROI and 4 - Layer
%
nps = cell(24,1);       % Number of pixels in each slice
tcpx = cell(24,1);      % Time constant - T1rho/T2* for each pixel
amppx = cell(24,1);     % Amplitude of exponential for each pixel
rsspx = cell(24,1);     % Sum of squared residuals for each pixel
%
% Loop through Compartments
%
for kc = 1:2
%
   ikc = kc-1;          % Compartment code (0 - lateral, 1 - medial)
   mskc = maskroi{kc};  % ROI masks for this compartment
   rslr = rsls{kc};     % Slices for this compartment
   nrslr = nrsls(kc);   % Number of slices in this compartment
%
% Loop through Bones
%
   for kb = 1:2
%
      ikb = kb-1;       % Bone code identifier (0 - femur, 1 - tibia)
      mskb = mskc{kb};  % ROI masks for this bone
      rslb = rslbs{kb}; % Slices for this bone
%
      maskb = masklay{kb};             % Layer mask for this bone
%
      rimgs = cell(nrslr,3,2);         % # of slices, # of ROIs, # of layers
%
% Loop through Slices
%
      np = zeros(3,2,nrslr);
%
      for ks = 1:nrslr
%
         slk = rslr(ks);               % Slice number
         idxl = rsl==slk;              % Index to layer slices in rsl/masklay
         idxb = rslb==slk;             % Index to bone slices in rslbs/maskroi
%
         rimgt = cell(ntime,3,2);      % # of spin lock/echo times, # of ROIs, # of layers
%
% Loop through Spin Lock/Echo Times
%
         for kt = 1:ntime
%
            rimg = squeeze(v(:,:,slk,kt));  % T1/T2 data for slice and spin lock/echo time
%
% Loop through Layers and ROIs
%
            for kl = 1:2               % 1 - superficial, 2 - deep
               maskl = maskb(:,kl,idxl);
               for kr = 1:3            % 1 - anterior/trochlea, 2 - central and 3 - posterior
                  mskr = mskb{kr};
                  msks = mskr(:,idxb);
                  mskl = msks&maskl;
                  rimgt{kt,kr,kl} = rimg(mskl)';
               end
            end
%
         end            % End of kt loop - Spin lock/echo times
%
         for kl = 1:2   % 1 - superficial, 2 - deep
            for kr = 1:3               % 1 - anterior/trochlea, 2 - central and 3 - posterior
               rimgs{ks,kr,kl} = cat(1,rimgt{:,kr,kl});    % Combine spin lock/echo times
               np(kr,kl,ks) = sum(mskl);
            end
         end
         
      end               % End of ks loop - slices loop
%
      clear rimg rimgt;
%
% Calculate T1rho/T2* for Each Layer
%
      for kl = 1:2      % 1 - superficial, 2 - deep
%
         ikl = 2-kl;    % 0 - deep, 1 -superficial
%
         for kr = 1:3   % 1 - anterior/trochlea, 2 - central and 3 - posterior
%
            idx = 12*kc+6*kb+3*kl+kr-21;    % Index to variables
            nps{idx} = np(kr,kl,:);         % Number of pixels in slice
%
            rimgl = cat(2,rimgs{:,kr,kl});
            npx(idx) = size(rimgl,2);
            r = rimgl(:);
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
            id(idx,:) = [ikc ikb kr-1 ikl];      % Compartment, bone, ROI, and layer IDs
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
         end            % End of kr loop - ROIs loop
%
      end               % End of kl loop - layers loop
   end                  % End of kb loop - bones loop
end                     % End of kc loop - compartments loop
%
return