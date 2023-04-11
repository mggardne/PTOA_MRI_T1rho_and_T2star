%#######################################################################
%
%            * CONTACT CHecK of PTOA Segmentations Program *
%
%          M-File which reads the segmentation CSV files to check if
%     the femur and tibia segmentations overlap.  If they do overlap,
%     the program measures the maximum amount of overlap.
%
%          The slices with overlap are plotted to Postscript files
%     contact_chk*.ps in the Results\ContactChecks folder.  The maximum
%     overlaps are also saved to the MS-Excel spreadsheet, 
%     contact_chk.xlsx, in the Results\ContactChecks folder.
%
%     NOTES:  1.  M-files lsect3.m, lsect4.m, lsect5.m, mxd2lins.m,
%             pts2lin.m, rd_prois.m, and rd_roi6.m must be in the
%             current directory or path.
%
%     15-Sep-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results','ContactChecks');    % Results directory
%
ifirst = true;          % First write to file
xlsnam = 'contact_chk.xlsx';           % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs = {'Subject' 'Result' 'Leg' 'Comprt' 'Slice', ...
        'MaxDist'};
%
psnam = fullfile(resdir,'contact_chk_');    % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Subject Directories
%
subjn = [3;7;12;17;18;19;20;21;24;30;31;32;34;35;38;46];
nsubj = size(subjn,1);
%
sdirs = int2str(subjn);
sdirs = [repmat('0',nsubj,1) sdirs];
sdirs = strrep(cellstr(sdirs),' ','0');
%
% Get Legs and Compartment Variables
%
legs = ['L'; 'R'];      % Left/Right
lnams = {'Left Leg'; 'Right Leg'};
%
cmprt  = {'Lateral'; 'Medial'};
%
% Initialize Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Result - 1 = T1rho and 2 = T2*
%   Index 3 - Leg - 1 = left and 2 = right
%   Index 4 - Compartment - 1 = lateral and 2 = medial
%
dmx = cell(nsubj,2,2,2);               % Maximum distances
sls = cell(nsubj,2,2,2);               % Slice numbers
dmxs = zeros(nsubj,2,2,2);             % Maximum distances for leg and compartment
%
% T1rho Segmentations
%
rdir = 'RHO';
ires = 1;               % T1rho
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory, Name and Number
%
   sdir = sdirs{ks};                   % Current subject directory
   subj = eval(sdir);                  % Subject number
%
   psnams = [psnam sdir];              % Add subject to PS file name
%
% T1rho Identifier
%
   psnamr = [psnams '_T1R_'];          % Add result type to PS file name
%
% Directory with Data Matlab MAT Files
%
   rdirk = fullfile(sdir,rdir);        % Directory with data
%
% Loop through Legs
%
   for kl = 1:2
%
% Get Leg
%
      leg = legs(kl);
      lnam = lnams{kl};
%
% Add Leg to PS File Name
%
      psnamf = [psnamr leg pstyp];     % Add leg to PS file name
      iplt1 = true;                    % First plot in new PS file
%
% Get CSV File
%
      csv_file = [sdir '_' leg '_ALL_' rdir '_*.csv'];
      dnam = dir(fullfile(rdirk,csv_file));
      dnam = {dnam.name}';
%
      if size(dnam,1)>1
        warning([' *** WARNING in contact_chkp: More than one', ...
                 ' segmentation file!  Using file:  ', dnam{1}]);
      end
%
% Read Segmentations
%
      brois = rd_prois(fullfile(rdirk,dnam{1}));
%
% Get Femur Data
%
      fc = brois(1).rois(1);           % Femur cartilage ROIs
      fsl = fc.slice;                  % All femur slices
      fdat = [fc.roi.data]';           % All femur data
%
% Tibia Compartment Data
%
      for m = 1:2 % 1 - lateral and 2 - medial
%
         tc = brois(2).rois(1);        % Tibia cartilage ROIs
%
         tcc = tc.roi(m);              % Tibia compartment
         tsl = tcc.imageno;            % Compartment slices
%
% Find Slices with Both Tibia and Femur Segmentations
%
         [~,idt,idf] = intersect(tsl,fsl);
%
% Loop through Compartment Slices
%
         ns = size(idt,1);
%
         for n = 1:ns
%
            fdats = fdat{idf(n)};
            tdats = tcc.data{idt(n)};
%
% Find Any Intersections Between the Tibia and Femur Segmentations
%
            [ipts,~,idcs] = lsect5(tdats,fdats);
%
            if ~isempty(ipts)
              figure;
              orient landscape;
%
              plot(fdats(:,1),fdats(:,2),'b.-','LineWidth', ...
                   1.5,'MarkerSize',12);
              hold on;
              plot(tdats(:,1),tdats(:,2),'g.-','LineWidth', ...
                   1.5,'MarkerSize',12);
              plot(ipts(:,1),ipts(:,2),'ro','LineWidth',1.5, ...
                   'MarkerSize',8);
              axis equal;
              set(gca,'XDir','reverse');
              set(gca,'YDir','reverse');
%
              title({['Subject ' sdir ' T1\rho']; ...
                     [cmprt{m} ' ' lnam ', Slice ' ...
                     int2str(tsl(idt(n)))]},'FontSize',16, ...
                    'FontWeight','bold');
%
              ni = size(ipts,1);
              dd = 0;
%
              for km = 1:ni/2
                 mi = 2*km;
                 idfx = fdats(:,1)>ipts(mi-1,1)&fdats(:,1)< ...
                              ipts(mi,1);
                 idtx = tdats(:,1)>ipts(mi-1,1)&tdats(:,1)< ...
                              ipts(mi,1);
                 nf = nnz(idfx);
                 nt = nnz(idtx);
                 df = 0;
                 dt = 0;
%
                 idcs(mi-1:mi,:) = sort(idcs(mi-1:mi,:));  % Make sure both lines are going in the same direction
%
                 if nf>0
                   idtl = idcs(mi-1,1):idcs(mi,1)+1;
                   txyz = tdats(idtl,:);
                   fxyz = fdats(idfx,:);
                   [df,idfmx,xyzf] = mxd2lins(txyz,fxyz);
                   if idfmx>0
                     plot([fxyz(idfmx,1); xyzf(:,1)], ...
                          [fxyz(idfmx,2); xyzf(:,2)],'r-', ...
                          'LineWidth',2);
                     text(xyzf(:,1),xyzf(:,2),sprintf('%.2f',df));
                   end
                 end
%
                 if nt>0
                   idfl = idcs(mi-1,2):idcs(mi,2)+1;
                   txyz = tdats(idtx,:);
                   fxyz = fdats(idfl,:);
                   [dt,idtmx,xyzt] = mxd2lins(fxyz,txyz);
                   if idtmx>0
                     plot([txyz(idtmx,1); xyzt(:,1)], ...
                          [txyz(idtmx,2); xyzt(:,2)],'m-', ...
                          'LineWidth',2);
                     text(xyzt(:,1),xyzt(:,2),sprintf('%.2f',dt));
                   end
                 end
%
% Get Slice Maximum Overlap
%
                 if df>dt&&df>dd
                   dd = df;
                 elseif dt>dd
                   dd = dt;
                 end
%
              end
%
              if iplt1
                print('-dpsc2','-r600','-fillpage',psnamf);
                iplt1 = false;
              else
                print('-dpsc2','-r600','-fillpage','-append', ...
                      psnamf);
              end
              close;
%
% Save Slice Maximum Overlap and Slice Number
%
              dmx{ks,ires,kl,m} = [dmx{ks,ires,kl,m}; dd];
              sls{ks,ires,kl,m} = [sls{ks,ires,kl, ...
                                        m}; fsl(idf(n))];
%
% Create and Write Table of Results
%
              dat = [subj ires-1 kl-1 m-1 fsl(idf(n)) dd];
              t = array2table(dat,'VariableNames',hdrs);
%
              if ifirst
                writetable(t,xlsnam,'WriteMode','replacefile');
                ifirst = false;
              else
                writetable(t,xlsnam,'WriteMode','append', ...
                           'WriteVariableNames',false);
              end
%
            end         % If femur and tibia segmentation intersect?
%
         end            % Slices loop - n
%
         if ~isempty(dmx{ks,ires,kl,m})
           dmxs(ks,ires,kl,m) = max(dmx{ks,ires,kl,m});
         end
%
      end               % End of tibia compartment loop - m
%
   end                  % End of leg loop - kl
%
end                     % End of subject loop - ks
%
% T2* Segmentations
%
rdir = 'T2S';
ires = 2;               % T2*
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory, Name and Number
%
   sdir = sdirs{ks};                   % Current subject directory
   subj = eval(sdir);                  % Subject number
%
   psnams = [psnam sdir];              % Add subject to PS file name
%
% T2* Identifier
%
   psnamr = [psnams '_T2S_'];          % Add result type to PS file name
%
% Directory with Data Matlab MAT Files
%
   rdirk = fullfile(sdir,rdir);        % Directory with data
%
% Loop through Legs
%
   for kl = 1:2
%
% Get Leg
%
      leg = legs(kl);
      lnam = lnams{kl};
%
% Add Leg to PS File Name
%
      psnamf = [psnamr leg pstyp];     % Add leg and load to PS file name
      iplt1 = true;                    % First plot in new PS file
%
% Get CSV File
%
      csv_file = [sdir '_' leg '_ALL_' rdir '_*.csv'];
      dnam = dir(fullfile(rdirk,csv_file));
      dnam = {dnam.name}';
%
      if size(dnam,1)>1
        warning([' *** WARNING in contact_chkp: More than one', ...
                 ' segmentation file!  Using file:  ', dnam{1}]);
      end
%
% Read Segmentations
%
      brois = rd_prois(fullfile(rdirk,dnam{1}));
%
% Get Femur Data
%
      fc = brois(1).rois(1);           % Femur cartilage ROIs
      fsl = fc.slice;                  % All femur slices
      fdat = [fc.roi.data]';           % All femur data
%
% Tibia Compartment Data
%
      for m = 1:2 % 1 - lateral and 2 - medial
%
         tc = brois(2).rois(1);        % Tibia cartilage ROIs
%
         tcc = tc.roi(m);              % Tibia compartment
         tsl = tcc.imageno;            % Compartment slices
%
% Find Slices with Both Tibia and Femur Segmentations
%
         [~,idt,idf] = intersect(tsl,fsl);
%
% Loop through Compartment Slices
%
         ns = size(idt,1);
%
         for n = 1:ns
%
            fdats = fdat{idf(n)};
            tdats = tcc.data{idt(n)};
%
% Find Any Intersections Between the Tibia and Femur Segmentations
%
            [ipts,~,idcs] = lsect5(tdats,fdats);
%
            if ~isempty(ipts)
              figure;
              orient landscape;
%
              plot(fdats(:,1),fdats(:,2),'b.-','LineWidth', ...
                   1.5,'MarkerSize',12);
              hold on;
              plot(tdats(:,1),tdats(:,2),'g.-','LineWidth', ...
                   1.5,'MarkerSize',12);
              plot(ipts(:,1),ipts(:,2),'ro','LineWidth',1.5, ...
                   'MarkerSize',8);
              axis equal;
              set(gca,'XDir','reverse');
              set(gca,'YDir','reverse');
%
              title({['Subject ' sdir ' T2*']; ...
                     [cmprt{m} ' ' lnam ', Slice ' ...
                     int2str(tsl(idt(n)))]},'FontSize',16, ...
                    'FontWeight','bold');
%
              ni = size(ipts,1);
              dd = 0;
%
              for km = 1:ni/2
                 mi = 2*km;
                 idfx = fdats(:,1)>ipts(mi-1,1)&fdats(:,1)< ...
                              ipts(mi,1);
                 idtx = tdats(:,1)>ipts(mi-1,1)&tdats(:,1)< ...
                              ipts(mi,1);
                 nf = nnz(idfx);
                 nt = nnz(idtx);
                 df = 0;
                 dt = 0;
%
                 idcs(mi-1:mi,:) = sort(idcs(mi-1:mi,:)); % Make sure both lines are going in the same direction
%
                 if nf>0
                   idtl = idcs(mi-1,1):idcs(mi,1)+1;
                   txyz = tdats(idtl,:);
                   fxyz = fdats(idfx,:);
                   [df,idfmx,xyzf] = mxd2lins(txyz,fxyz);
                   if idfmx>0
                     plot([fxyz(idfmx,1); xyzf(:,1)], ...
                          [fxyz(idfmx,2); xyzf(:,2)],'r-', ...
                          'LineWidth',2);
                     text(xyzf(:,1),xyzf(:,2),sprintf('%.2f',df));
                   end
                 end
%
                 if nt>0
                   idfl = idcs(mi-1,2):idcs(mi,2)+1;
                   txyz = tdats(idtx,:);
                   fxyz = fdats(idfl,:);
                   [dt,idtmx,xyzt] = mxd2lins(fxyz,txyz);
                   if idtmx>0
                     plot([txyz(idtmx,1); xyzt(:,1)], ...
                          [txyz(idtmx,2); xyzt(:,2)],'m-', ...
                          'LineWidth',2);
                     text(xyzt(:,1),xyzt(:,2),sprintf('%.2f',dt));
                   end
                 end
%
% Get Slice Maximum Overlap
%
                 if df>dt&&df>dd
                   dd = df;
                 elseif dt>dd
                   dd = dt;
                 end
%
              end
%
              if iplt1
                print('-dpsc2','-r600','-fillpage',psnamf);
                iplt1 = false;
              else
                print('-dpsc2','-r600','-fillpage','-append', ...
                      psnamf);
              end
              close;
%
% Save Slice Maximum Overlap and Slice Number
%
              dmx{ks,ires,kl,m} = [dmx{ks,ires,kl,m}; dd];
              sls{ks,ires,kl,m} = [sls{ks,ires,kl, ...
                                        m}; fsl(idf(n))];
%
% Create and Write Table of Results
%
              dat = [subj ires-1 kl-1 m-1 fsl(idf(n)) dd];
              t = array2table(dat,'VariableNames',hdrs);
%
              if ifirst
                writetable(t,xlsnam,'WriteMode','replacefile');
                ifirst = false;
              else
                writetable(t,xlsnam,'WriteMode','append', ...
                           'WriteVariableNames',false);
              end
%
            end         % If femur and tibia segmentation intersect?
%
         end            % Slices loop - n
%
         if ~isempty(dmx{ks,ires,kl,m})
           dmxs(ks,ires,kl,m) = max(dmx{ks,ires,kl,m});
         end
%
      end               % End of tibia compartment loop - m
%
   end                  % End of leg loop - kl
%
end                     % End of subject loop - ks
%
% Get Overlaps by Analysis Type
%
dmx1 = squeeze(dmx(:,:,1,:,:,:));      % T1rho
dmx2 = squeeze(dmx(:,:,2,:,:,:));      % T2*
%
dmx1 = cell2mat({dmx1{:}}');
dmx2 = cell2mat({dmx2{:}}');
%
dmx1 = dmx1(dmx1>0);
dmx2 = dmx2(dmx2>0);
%
% Overlap Descriptive Statistics
%
n1 = size(dmx1,1);
n2 = size(dmx2,1);
%
n21 = nnz(dmx1>2);      % Overlaps > 2 pixels
n22 = nnz(dmx2>2);      % Overlaps > 2 pixels
%
p21 = 100*n21/n1;
p22 = 100*n22/n2;
%
mn1 = mean(dmx1);
mx1 = max(dmx1);
sd1 = std(dmx1);
mn2 = mean(dmx2);
mx2 = max(dmx2);
sd2 = std(dmx2);
%
% Overlap Histograms
%
figure;
orient landscape;
histogram(dmx1,'BinWidth',0.05);
xlabel('Overlap (pixels)','FontSize',12,'FontWeight','bold');
ylabel('Frequency','FontSize',12,'FontWeight','bold');
title('T1\rho Overlaps','FontSize',16,'FontWeight','bold');
%
axlim = axis;
hold on;
% plot([2 2],axlim(3:4),'k:','LineWidth',0.5);
ptxt{3} = sprintf('SD = %.2f pixels',sd1);
ptxt{2} = sprintf('Maximum = %.2f pixels',mx1);
ptxt{1} = sprintf('Mean = %.2f pixels',mn1);
text(0.55,sum(axlim(3:4))/2,ptxt,'FontSize',12,'FontWeight','bold');
ptxt = sprintf('%.2f%% of overlaps > 2 pixels',p21);
ptxt = {[int2str(n21) ' overlaps > 2 pixels']; ptxt};
text(0.8,sum(axlim(3:4))/4,ptxt,'FontSize',12,'FontWeight','bold');
%
psnamh = fullfile(resdir,'overlap_hist.ps');
print('-dpsc2','-r600','-fillpage',psnamh);
%
figure;
orient landscape;
histogram(dmx2,'BinWidth',0.1);
xlabel('Overlap (pixels)','FontSize',12,'FontWeight','bold');
ylabel('Frequency','FontSize',12,'FontWeight','bold');
title('T2* Overlaps','FontSize',16,'FontWeight','bold');
%
axlim = axis;
hold on;
% plot([2 2],axlim(3:4),'k:','LineWidth',0.5);
ptxt{3} = sprintf('SD = %.2f pixels',sd2);
ptxt{2} = sprintf('Maximum = %.2f pixels',mx2);
ptxt{1} = sprintf('Mean = %.2f pixels',mn2);
text(0.9,sum(axlim(3:4))/2,ptxt,'FontSize',12,'FontWeight','bold');
ptxt = sprintf('%.2f%% of overlaps > 2 pixels',p22);
ptxt = {[int2str(n22) ' overlaps > 2 pixels']; ptxt};
text(1.4,sum(axlim(3:4))/4,ptxt,'FontSize',12,'FontWeight','bold');
print('-dpsc2','-r600','-append','-fillpage',psnamh);
%
return