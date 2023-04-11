function [xyzo,xyzax,xyzlnx,r,s,sc,sy,ssep,hf] = f_cs_14(datlb, ...
                                             datmb,datfs,datcp,iplt,hf)
%F_CS_14   Uses MRI femur (F) data cell arrays to calculate a femur
%          based coordinate system (CS).
%
%          [XYZO,XYZAX] = F_CS_14(DATLB,DATMB,DATFS,DATCP)
%          given the MRI femur data cell arrays for the lateral condyle,
%          DATLB, the medial condyle, DATMB, the femoral shaft, DATFS
%          and the center point of the notch, DATCP, calculates a femur
%          coordinate system based on a cylinder fit of the condyle bone
%          surfaces, the centroid of the femur shaft from just above
%          the patella and the center point of the notch.  The function
%          returns the origin of the femur coordinate system, in XYZO,
%          and unit X, Y and Z vectors for the X, Y, and Z axes in the
%          rows of matrix XYZAX (rotation matrix).
%
%          [XYZO,XYZAX,R] = F_CS_14(DATLB,DATMB,DATFS,DATCP) returns
%          the radius of the fitted cylinder, R.
%
%          NOTES:  1.  The Matlab files cyl_fit.m, cyl_plt.m,
%                  dist2cyl.m, plt_datsl.m, pts2lin.m, pt2line.m, and
%                  tri_area.m must be in the current directory or path.
%
%                  2.  Note the femur coordinate system is defined as
%                  the positive X-axis as either anterior (left knees)
%                  or posterior (right knees), positive Y-axis as
%                  lateral and positive Z-axis as superior.
%
%                  3.  Based on the femur coordinate system work of 
%                  Daniel Sturnick.
%
%          29-Jul-2014 * Mack Gardner-Morse
%
%          29-Nov-2022 * Mack Gardner-Morse * Uses the new cyl_fit.m
%          function which uses the Matlab nonlinear least squares
%          function lsqnonlin for fitting the cylinder.
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  error(' *** ERROR in F_CS_14:  Four input variables are required!');
end
%
if (nargin<5)||isempty(iplt)
  iplt = false;
end
%
if (nargin<6)
    if iplt
      hf = figure;
      orient landscape;
      hold on;
      view(3);
    end
else
  figure(hf);
  orient landscape;
  hold on;
  view(3);
end
%
% Read Femoral Shaft Data and Get Scaling
%
xyzfs = datfs{1};
nptfs = size(xyzfs,1);
%
xyzm = mean(xyzfs);
xyzp = [xyzfs; xyzm];
%
if iplt
  plot3(xyzfs(:,1),xyzfs(:,2),xyzfs(:,3),'k.-','MarkerSize',6, ...
        'LineWidth',1);
end
%
s = max(xyzfs)-min(xyzfs);             % Scaling in MRI CS
%
% Get and Plot Centroid
%
tri = [repmat(nptfs+1,nptfs-1,1) (1:nptfs-1)' (2:nptfs)'];
tri = [tri; [nptfs+1 nptfs 1]];
%
[at,cgt] = tri_area(xyzp(:,1),xyzp(:,2),xyzp(:,3),tri);
xyzc = sum(repmat(at,1,3).*cgt)./sum(at);
%
if iplt
  plot3(xyzc(1),xyzc(2),xyzc(3),'bs','LineWidth',2,'MarkerSize',6);
end
%
% Get Notch Coordinates
%
xyz0 = datcp{1};
if iplt
  plot3(xyz0(1),xyz0(2),xyz0(3),'gs','Color',[0 0.5 0], ...
        'LineWidth',2,'MarkerSize',6);
end
%
% Plot a Line from Notch to Centroid of Slice 1 (Z-Axis Direction)
%
xyzl2 = [xyz0; xyzc];
if iplt
  plot3(xyzl2(:,1),xyzl2(:,2),xyzl2(:,3),'g-','Color',[0 0.5 0], ...
        'LineWidth',2);
end
%
% Add Labels to Figure
%
if iplt
  xlabel('X (mm)','FontSize',12,'FontWeight','bold');
  ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
  zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
%   title(title_txt,'FontSize',16,'FontWeight','bold','interpreter','none');
  axis equal;
end   
%
% Get Condyle Sagittal Data
%
xyzlat = cell2mat(datlb);
xyzmed = cell2mat(datmb);
%
% [nslice,nsl,isl] = sl_info(datlb);
%
if iplt
  plt_datsl(datlb,'b.-',0.5,[0 0 0.7]);
  plt_datsl(datmb,'g.-',0.5,[0 0.5 0]);
end
%
% Fit Cylinder to Slice Data
%
ri = 15;                % Near minimum radius
%ri = 25;                % For Knee ID #'s 8,30,62,71
%ri = 32                 % For Knee ID #'s 29,40  
xyz1i = mean(xyzlat);
xyz2i = mean(xyzmed);
[r,xyz1,xyz2,ssep] = cyl_fit(ri,xyz1i,xyz2i,[xyzlat; xyzmed]);
%
xyzlnx = [xyz1; xyz2];
%
if iplt
  hc = cyl_plt(r,xyz1,xyz2);
  set(hc,'EdgeColor',[0.5 0.5 0.5]);
%
  plot3(xyzlnx(:,1),xyzlnx(:,2),xyzlnx(:,3),'ko-','LineWidth',2, ...
        'Color',[0.5 0.5 0.5]);
end
%
% Get Coordinate System
%
yax = -diff(xyzlnx);
sy = norm(yax);         % Width of medial and lateral points along Y-axis
yax = yax./sy;
zax = diff(xyzl2);
zax = zax./norm(zax);
% zax = sign(zax(3))*rvec;
xax = cross(yax,zax);
xax = xax./norm(xax);
zax = cross(xax,yax);
zax = zax./norm(zax);
%
% Get Axis Origin on Cylinder Axis Line Nearest Notch Point
% (pt2line.m)
%
xyzo = pt2line(xyz1,yax,xyz0);
%
xyzax = [xax; yax; zax]';              % Rotation matrix
%
% Plot Coordinate System Origins
%
if iplt
  xc = repmat(xyzo(1),1,3);
  yc = repmat(xyzo(2),1,3);
  zc = repmat(xyzo(3),1,3);
%
  plot3(xc(1),yc(1),zc(1),'ro','MarkerSize',9,'LineWidth',2);
  h = quiver3(xc,yc,zc,xyzax(1,:),xyzax(2,:),xyzax(3,:),50,'c'); 
  set(h,'LineWidth',3);
%
  axis equal;
end
%
% Femoral Shaft Scaling
%
xyzfst = xyzfs-repmat(xyzo,nptfs,1);
xyzfst = xyzfst*xyzax';
%
sc = max(xyzfst)-min(xyzfst);          % Scaling in femur CS
%
return