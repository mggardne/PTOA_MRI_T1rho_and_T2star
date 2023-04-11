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
