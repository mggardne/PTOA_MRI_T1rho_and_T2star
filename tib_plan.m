function [lcap,lcav,lpcp,lpcv,mcap,mcav,mpcp,mpcv] = tib_plan(t3, ...
                                                 idl,idm,iplt,ltxt,fsav)
%TIB_PLAN  Uses tibia bone sagittal plane segmentations to calculate
%          planes dividing the tibia into regions of interest (ROIs).
%
%          [LCAP,LCAV,LPCP,LPCV,MCAP,MCAV,MPCP,MPCV] = TIB_PLAN(T3,IDL,
%          IDM) Given a sagittal knee segmentation cell array, T3,
%          index to lateral compartment slices, IDL, and index to
%          medial slices, IDM, calculates the lateral central-anterior
%          dividing plane defined by a point, LCAP, and a normal vector,
%          LMV, the lateral posterior-central dividing plane defined by
%          a point, LPCP, and a normal vector, LPCV, the medial central-
%          anterior dividing defined by a point, MCAP, and a normal
%          vector, MPCP, and a medial posterior-central plane defined
%          by a point, MPCP, and a normal vector, MPCV.
%
%          [...] = TIB_PLAN(F3,IDL,IDM,IPLT,LTXT,FSAV) given a logical
%          true, IPLT, creates plots of the tibial dividing planes.
%          LTXT is a string with the name of the leg to be used in the
%          plot title.  If the file name, FSAV, is provided, the plots
%          are appended to the file, FSAV.
%
%          NOTES:  1.  The tibia bone sagittal plane segmentations
%                  should be in the tibial coordinate system.
%
%                  2.  The Matlab files plnorm.m and plt_datsl.m must
%                  be in the current directory or path.
%
%                  3.  The plot is appended to the file, FSAV.  If the
%                  does not exist, Matlab will create the file.
%
%          17-Jan-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error([' *** ERROR in TIB_PLAN:  Sagittal tibial segmentations', ...
         ' and indices to lateral and medial compartments are', ...
         ' required inputs!']);
end
%
if (nargin<4)||isempty(iplt)
  iplt = false;
else
  iplt = true;
end
%
if (nargin<5)||isempty(ltxt)
  ltxt = [];
end
%
if (nargin<6)||isempty(fsav)
  iprt = false;
else
  iprt = true;
end
%
% X-axis is Normal Vector for All Planes
%
[lcav,lpcv,mcav,mpcv] = deal([1 0 0]);
%
% Get Anterior-Posterior Length of the Lateral and Medial Compartments
%
xyzl = cell2mat(t3(idl));
ctrl = mean(xyzl);      % Center of lateral points
minlatx = min(xyzl(:,1));
maxlatx = max(xyzl(:,1));
dlat = maxlatx-minlatx;
%
xyzm = cell2mat(t3(idm));
ctrm = mean(xyzm);      % Center of medial points
minmedx = min(xyzm(:,1));
maxmedx = max(xyzm(:,1));
dmed = maxmedx-minmedx;
%
% Get Points in the Dividing Planes
%
dlat3 = dlat/3;         % Divide distance into thirds
lcap = [minlatx+2*dlat3 ctrl(:,2:3)];
lpcp = [minlatx+dlat3 ctrl(:,2:3)];
%
dmed3 = dmed/3;         % Divide distance into thirds
mcap = [minmedx+2*dmed3 ctrm(:,2:3)];
mpcp = [minmedx+dmed3 ctrm(:,2:3)];
%
if iplt
  figure;
  orient landscape;
  plt_datsl(t3(idl),'g.-',0.5,[0 0.7 0]);   % Lateral
  hold on;
  plt_datsl(t3(idm),'b.-',0.5,[0 0 0.7]);   % Medial
  axis equal;
  view(20,32);
  grid on;
  xlabel('X (mm)','FontSize',12,'FontWeight','bold');
  ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
  zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
  title({'Dividing Planes in Tibia CS'; ltxt},'FontSize',16, ...
        'FontWeight','bold');
  axis equal;
  %
  axlim = axis;
  %
  pxyz1 = [lcap(1) 0 axlim(5); lcap(1) 0 axlim(6); ...
          lcap(1) axlim(4) axlim(6); lcap(1) axlim(4) axlim(5); ...
          lcap(1) 0 axlim(5)];
  patch(pxyz1(:,1),pxyz1(:,2),pxyz1(:,3),pxyz1(:,3),'FaceColor', ...
        [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
  %
  pxyz2 = [lpcp(1) 0 axlim(5); lpcp(1) 0 axlim(6); ...
          lpcp(1) axlim(4) axlim(6); lpcp(1) axlim(4) axlim(5); ...
          lpcp(1) 0 axlim(5)];
  patch(pxyz2(:,1),pxyz2(:,2),pxyz2(:,3),pxyz2(:,3),'FaceColor', ...
        [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
  %
  pxyz3 = [mcap(1) axlim(3) axlim(5); mcap(1) axlim(3) axlim(6); ...
          mcap(1) 0 axlim(6); mcap(1) 0 axlim(5); ...
          mcap(1) axlim(3) axlim(5)];
  patch(pxyz3(:,1),pxyz3(:,2),pxyz3(:,3),pxyz3(:,3),'FaceColor', ...
        [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
  %
  pxyz4 = [mpcp(1) axlim(3) axlim(5); mpcp(1) axlim(3) axlim(6); ...
          mpcp(1) 0 axlim(6); mpcp(1) 0 axlim(5); ...
          mpcp(1) axlim(3) axlim(5)];
  patch(pxyz4(:,1),pxyz4(:,2),pxyz4(:,3),pxyz4(:,3),'FaceColor', ...
        [0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33);
%
  if iprt
    print('-dpsc2','-r600','-fillpage',fsav);
    view(-90,90);
    print('-dpsc2','-r600','-fillpage','-append',fsav);
  end
%
end
%
return