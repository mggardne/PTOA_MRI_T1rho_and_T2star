%#######################################################################
%
%                      * COMBine LAYers Program *
%
%          M-file which reads the PTOA qMRI results spreadsheet and 
%     combines the results for the two layers (superficial and deep).
%     The original and combined results are output to a spreadsheet
%     with all three results (0 - superficial, 1 - deep, and 2 - both
%     layers.  The new spreadsheet name starts with the old spreadsheet
%     name and appends "_both" to the name.
%
%     08-Jul-2026 * Mack Gardner-Morse
%

%#######################################################################
%
% Get PTOA qMRI Results MS-Excel Spreadsheet File
%
[snam,pnam] = uigetfile({'mri_*.xlsx', ...
                'qMRI spreadsheet'}, ...
                'Please Select a Spreadsheet File');
%
if isequal(snam,0)
  return;
end
%
fnam = fullfile(pnam,snam);
%
% Read Spreadsheet 
%
t = readtable(fnam);
%
% Index Variables
%
id_vars = {'Subject','Result','Leg','Bone','Comprt','ROI','Layer'};
nid = size(id_vars,2);
res_vars = {'ValidPix','Mean'};
nres = 2;
vars = [id_vars res_vars];
nvars = size(vars,2);
%
% Get Size of Table
%
nrows = size(t,1);
nlays = 2;
nboth = nrows/nlays;
%
% Sort Rows of Table
%
t = sortrows(t,id_vars);
%
% Get ID and Results Columns
%
tl = reshape(t{:,vars}',nvars,nlays,nboth);
%
% Check for Missing Layers
%
idl = 1:nid-1;
id = find(any(tl(idl,1,:)~=tl(idl,2,:)));
if ~isempty(id)
  error(' *** ERROR in comb_lay:  Missing layer!');
end
%
% Combine Layers
%
tl0 = squeeze(tl(end-1:end,1,:));
tl1 = squeeze(tl(end-1:end,2,:));
n0 = tl0(1,:);
m0 = tl0(2,:);
n1 = tl1(1,:);
m1 = tl1(2,:);
n2 = n0+n1;
idv = n0&n1;
res2 = (n0.*m0+n1.*m1)./n2;
res2(~idv) = 0;
%
% Create Combined Layers Table
%
t2 = t(1:2:end,:);
t2.ValidPix = n2';
t2.Mean = res2';
t2.Layer = repmat(2,nboth,1);
%
t2{:,{'Min','Max','SD','COV'}} = zeros(nboth,4);
%
t2 = [t; t2];
t2 = sortrows(t2,id_vars);
%
% Write to a New Spreadsheet
%
idot = strfind(fnam,'.xls');
fnam2 = [fnam(1:idot-1) '_both' fnam(idot:end)];
writetable(t2,fnam2);
%
return