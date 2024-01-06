function [pulpData, pulpDataFiltered] = importPulpCharacterizationData(characterizationDataDir, characterizationDataFileName, MyPackingFilterSettings)

rawData = dlmread(horzcat(characterizationDataDir,filesep,characterizationDataFileName));

% Data is structured as follows:
% 2. Lp [mm]
% 3. Lc [mm]
% 4. Width [um]
% 5. Wall thickness [um]
% 6. Curl [%]

pulpData.Lp = rawData(:,2)  / 1e3;
pulpData.Lc = rawData(:,3) / 1e3;
pulpData.Width = rawData(:,4) / 1e6;
pulpData.WallThickness = rawData(:,5) / 1e6;
pulpData.Curl = rawData(:,6) ./ 100.0;


% Apply filters
filteredData = rawData;
filteredData(filteredData(:,2)/ 1e3 < MyPackingFilterSettings.Lc(1) , :) = [];
filteredData(filteredData(:,2)/ 1e3 > MyPackingFilterSettings.Lc(2) , :) = [];

filteredData(filteredData(:,6)/100.0 < MyPackingFilterSettings.Curl(1) , :) = [];
filteredData(filteredData(:,6)/100.0 > MyPackingFilterSettings.Curl(2) , :) = [];

filteredData(:,4) = filteredData(:,4) * MyPackingFilterSettings.Width(3);

filteredData(:,5) = filteredData(:,5) * MyPackingFilterSettings.WallThickness(3);

pulpDataFiltered.Lp = filteredData(:,2)  / 1e3;
pulpDataFiltered.Lc = filteredData(:,3) / 1e3;
pulpDataFiltered.Width = filteredData(:,4) / 1e6;
pulpDataFiltered.WallThickness = filteredData(:,5) / 1e6;
pulpDataFiltered.Curl = filteredData(:,6) ./ 100.0;


disp(['Filtering reduced the observations from ' num2str(size(rawData,1)) ' to ' num2str(size(filteredData,1))])