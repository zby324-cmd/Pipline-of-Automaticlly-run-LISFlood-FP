function [sh,eh] = FERCallerS(csvHydrograph,strBasin,dirCPM)
curFile = mfilename('fullpath');
[curDir,~,~]=fileparts(curFile);
[progDir,~,~]=fileparts(curDir);
addpath([progDir,'/CREST_Prep']);
addpath([progDir,'/MEX']);
% addpath([progDir,'/CREST_3.0/common']);
addpath([progDir,'/geomorphology/common']);
addpath([progDir,'/Hydrograph_separation']);
addpath([progDir,'/ThirdPartyM/nansuite']);
addpath([progDir,'/ThirdPartyM/statbx42']);
addpath([progDir,'/ThirdPartyM/xticklabel_rotate']);
addpath([progDir,'/ThirdPartyM/TDF']);
addpath([progDir,'/ThirdPartyM/rdir']);
GDALLoad();
%% read flow and precip
% [data,txt]=xlsread(csvHydrograph);
% dates=txt(2:end,1);

% read csv file
dataTable = readtable(csvHydrograph);

% extract time and data 
dates = dataTable{:, 1}; 
data = table2array(dataTable(:, 2:end));

dn=datenum(dates);
clear txt dates
QFilled=data(:,end);
clear data
QFilled(QFilled<0)=NaN;
%% Adjust QFilled size to 24*N  - Lzq25
lenQ = length(QFilled);
remainder = mod(lenQ, 24);
if remainder ~= 0
   fprintf('origin data records length are %d, not 24*N. \n', lenQ);
   % two ways:fill or cut
   % way 1 : cut
   %QFilled = QFilled(1:end-remainder);
   %dn = dn(1:end-remainder);
   %fprintf('cut the records, new length is  %d。\n', length(QFilled));
        
   % way2: fill NaN
   padding = 24 - remainder;
   QFilled = [QFilled; NaN(padding, 1)];
   dn = [dn; NaN(padding, 1)];
   fprintf('Filled records, new length is  %d. \n', length(QFilled));
end
%% read basin info
[raster,geoTrans,proj]=ReadRaster(strBasin);
bGCS=IsGeographic(proj,geoTrans);
area=RasterArea(raster,geoTrans,bGCS);
%% extraction
Recession_Analysis2(QFilled,[],area,dirCPM);
Baseflow_Separation3(dirCPM);
[sh,eh]=Events_Formation4_1(dirCPM,80);
end
