
path(path,'/shared/stormcenter/Linzq25/B002_Project_DOE/Z999_Tool')
path(path,'/shared/stormcenter/Linzq25/B002_Project_DOE/Z999_Tool/CREST_Prep')
path(path,'/shared/stormcenter/Linzq25/B002_Project_DOE/Z999_Tool/FloodEvent')

clear;clc;

savepath = '/shared/stormcenter/Linzq25/B002_Project_DOE/A001_USGS_sta_TVA/Transfer_hourly_Fix_format/';
csvFiles = dir(fullfile(savepath, '*.csv'));
USGS_fileNames = cell({csvFiles.name});

strBasin = '/shared/stormcenter/Linzq25/B002_Project_DOE/A002_basic/basinMask.tif';

Oufolder= '/shared/stormcenter/Linzq25/B002_Project_DOE/A003_extract_Flood_event_fix_format/';
if ~exist(Oufolder,'dir')
    mkdir(Oufolder);
end

for i = 231 : length(USGS_fileNames)

	%input xlsx file
	csvHydrograph = [savepath,USGS_fileNames{i}];

	% output
	dirCPM = [Oufolder,USGS_fileNames{i}(1:end-6),'/'];

	% mkdir
	if ~exist(dirCPM, 'dir')
		mkdir(dirCPM);
	end
	cd('/shared/stormcenter/Linzq25/B002_Project_DOE/Z999_Tool/FloodEvent')
	[sh,eh] = FERCallerS(csvHydrograph,strBasin,dirCPM);

	save([dirCPM,'Start_and_End_Time.mat'],'sh','eh')
	disp([num2str(i),' ',USGS_fileNames{i}(1:end-5),' is done!'])
	clear dirCPM csvHydrograph sh eh

end

