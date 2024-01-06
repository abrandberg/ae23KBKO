% Post processing of MyPacking output files.
%
% ABOUT:
%       - Takes generated MyPacking.exe outputs and generates a summary of the
%         fiber shapes.
%
% created:  06-01-2024
% author:   August Brandberg august.brandberg@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 01. Initialize
clear
clc
close all
format compact
format long

legendSettings = {'Characterization file','Char. f. after filtering','MyPacking','location','northeast'};
pulpDataHistSettings = {'facecolor','b', 'facealpha', 0.25};
pulpDataFilteredHistSettings = {'facecolor','r', 'facealpha', 0.25};
myPackingHistSettings = {'facecolor','g', 'facealpha', 0.25};

% 02. Inputs
workDir = '';
% Name of the FOLDER containing MyPacking files are located.

networkName = '';
% Name of the files, without extension.

truncationXDIM = 0.003*[-inf inf]; % If inf, look at all data.
truncationYDIM = 0.003*[-inf inf]; % If inf, look at all data.
% Sets limits of what to consider. Can be useful when evaluating e.g., mass of the network, as
% density might drop off close to the outer edges of the network (depends on the exact version
% of MyPacking.exe being used). Set to [-inf inf] to consider everything.


characterizationDataDir = workDir;
% Name of the FOLDER containing pulp characterization data used as input to MyPacking.exe
characterizationDataFileName = '';
% Name of the pulp characterization file, with extension.

% Value 1: Lower limit
% Value 2: Upper limit
% Value 3: Linear scaling factor
MyPackingFilterSettings.Lc            = [0.2d-3  inf    0   ];
MyPackingFilterSettings.Curl          = [0.0     2.09   0   ];
MyPackingFilterSettings.Width         = [-inf    inf    0.78];
MyPackingFilterSettings.WallThickness = [-inf    inf    0.53];


% 03. Import data
[nodalData,realSetData,elementData] = importNetworks(workDir,networkName);
% MyPacking data import

[pulpData, pulpDataFiltered] = importPulpCharacterizationData(characterizationDataDir, characterizationDataFileName, MyPackingFilterSettings);
% Characterization data import

% 04. Process data and generate effective properties such as fiber length.
networkStructure = tabulateNetworks(nodalData,elementData,realSetData,[truncationXDIM ; truncationYDIM],networkName,workDir);

% 05. Plot results
% Fig A: Polar plot of in-plane fiber segment orientations
figure
%polarhistogram(networkStructure(selIdx).angYX,360,'displaystyle','bar','Normalization','pdf');
rose(networkStructure.angYX,360);


% Fig B: Check anisotropy through the depth of the sheet
zDim = [min(networkStructure.meanPos(:,3)) max(networkStructure.meanPos(:,3))];
colorTemp = winter(10);
figure
for tLoop = 0:9
    selIdx =   networkStructure.meanPos(:,3) > (zDim(1) + 0.1*tLoop*(zDim(2)-zDim(1)))  ...
             & networkStructure.meanPos(:,3) < (zDim(1) + (0.1+0.1*tLoop)*(zDim(2)-zDim(1)));
    %polarhistogram(networkStructure(1).angYX(selIdx),360,'displaystyle','stairs','Normalization','pdf','edgeColor',colorTemp(tLoop,:));
    rose(networkStructure.angYX(selIdx),360);
    hold on

    drawnow;
    disp(['Number of segments in bin: ' num2str(tLoop) ' : '  num2str(sum(selIdx))])
end

% Fig C: Fiber wall thickness
% The types are: 1. Solid rectangle (corresponding to closed lumen)
%                2. Hollow rectangle (open lumen)
% Fibers with closed lumen have "0.0" indicated as wall thickness.
% Fibers are always oriented so that width >= height.
% Therefore, we can take the thickness of fiber wall in the solid rectangles as
% height * 0.5

selIdxSolid = networkStructure.realSetData(:,5) == 0.0;
selIdxHollow = not(selIdxSolid);
disp(['Solid cross-sections : ' num2str(sum(selIdxSolid)) ' , Hollow cross-sections : ' num2str(sum(selIdxHollow))])

fiberThickness(selIdxSolid) = 0.5*networkStructure.realSetData(selIdxSolid,4);
fiberThickness(selIdxHollow) = networkStructure.realSetData(selIdxHollow,5);

figure
figureBinVector = linspace(min(fiberThickness),max(fiberThickness),25);
hist(pulpData.WallThickness,figureBinVector,1,pulpDataHistSettings{:})
hold on
hist(pulpDataFiltered.WallThickness,figureBinVector,1,pulpDataFilteredHistSettings{:})
hist(fiberThickness(selIdxHollow),figureBinVector,1,myPackingHistSettings{:})
hist(fiberThickness(selIdxSolid),figureBinVector,1,'facecolor','m','facealpha',0.25)
xlabel('Fiber wall thickness [m]')
ylabel('Relative frequency [-]')
legend('Characterization file','Char. f. after filtering','MyPacking (hollows)', 'MyPacking (solids)', 'location','northeast')


% Fig D: Fiber width
figure
figureBinVector = linspace(min(networkStructure.realSetData(:,3)),max(networkStructure.realSetData(:,3)),25);
hist(pulpData.Width,figureBinVector,1,pulpDataHistSettings{:})
hold on
hist(pulpDataFiltered.Width,figureBinVector,1,pulpDataFilteredHistSettings{:})
hist(networkStructure.realSetData(:,3),figureBinVector,1,myPackingHistSettings{:})
xlabel('Fiber width [m]')
ylabel('Relative frequency [-]')
legend(legendSettings{:})

% Fig E: Fiber height
figure
figureBinVector = linspace(min(networkStructure.realSetData(:,4)),max(networkStructure.realSetData(:,4)),25);
hist(networkStructure.realSetData(:,4),figureBinVector,1,myPackingHistSettings{:})
xlabel('Fiber height [m]')
ylabel('Relative frequency [-]')
%legend(legendSettings{:})

% Fig F: Fiber projected length
figure
figureBinVector = linspace(min(networkStructure.fiberLengthProjected),max(networkStructure.fiberLengthProjected),25);
hist(pulpData.Lp,figureBinVector,1,pulpDataHistSettings{:})
hold on
hist(pulpDataFiltered.Lp,figureBinVector,1,pulpDataFilteredHistSettings{:})
hist(networkStructure.fiberLengthProjected,figureBinVector,1,myPackingHistSettings{:})
xlabel('Fiber projected length [m]')
ylabel('Relative frequency [-]')
legend(legendSettings{:})

% Fig G: Fiber length
figure
figureBinVector = linspace(min(networkStructure.fiberLength),max(networkStructure.fiberLength),25);
hist(pulpData.Lc,figureBinVector,1, pulpDataHistSettings{:})
hold on
hist(pulpDataFiltered.Lc,figureBinVector,1,pulpDataFilteredHistSettings{:})
hist(networkStructure.fiberLength,figureBinVector,1,myPackingHistSettings{:})
xlabel('Fiber length [m]')
ylabel('Relative frequency [-]')
legend(legendSettings{:})

% Fig H: Fiber curl
figure
figureBinVector = linspace(min(networkStructure.fiberCurl),max(networkStructure.fiberCurl),25);
hist(pulpData.Curl,figureBinVector,1, pulpDataHistSettings{:})
hold on
hist(pulpDataFiltered.Curl,figureBinVector,1,pulpDataFilteredHistSettings{:})
hist(networkStructure.fiberCurl,figureBinVector,1,myPackingHistSettings{:})
xlabel('Fiber curl [-]')
ylabel('Relative frequency [-]')
legend(legendSettings{:})




