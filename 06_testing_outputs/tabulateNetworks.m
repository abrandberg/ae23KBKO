function networkStructure = tabulateNetworks(nodalData,elementData,realSetData,probeSpace,networkName,networkDir)
%function
%tabulateNetworks(nodalData,elementData,realSetData,probeSpace,networkName,networkDir)
%uses the fibnet data to calculate various things. In this version, the
%characterization of interest is the fiber orientation.
%
% INPUT:    nodalData   : All nodal coordinates {IDX X Y Z}
%           elementData : All element connectivity {IDX START_NODE END_NODE MID_NODE REAL_IDX MAT_IDX}
%           realSetData : Real values {IDX TYPE WIDTH HEIGHT WALL_TKN}
%                         TYPE = 1 : Solid rectangle
%                         TYPE = 2 : Hollow rectangle
%           probeSpace  : Rectangle specifying region of interest ([xMin xMax ; yMin yMax])
%           networkDir  : The directory to import from
%           networkName : A string matching the network to import
%
% OUTPUT:   networkStructure : Calculated data, in this case fiber segment
%                              orientation. Also contains some fields for
%                              filtering purposes.
%
% ABOUT:
%
% TO DO:
% - This script does not bother itself with partial fiber segments (where 1
%   part is inside and 1 part is outside the region of interest). I have
%   not checked but I do not think it matters. /AB, 09-07-2019
%
% created: 09-07-2019
% author: August Brandberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xDim = probeSpace(1,:); % Boundaries for the region of interest
yDim = probeSpace(2,:);
zDim = [min(nodalData(:,4)) max(nodalData(:,4))];


% Pre-allocation
newFiberMarker = diff(elementData(:,5));
newFiberMarker = [1 ; newFiberMarker];
newFiberIndex = find(newFiberMarker);
angYX = nan(size(elementData(:,1)));%[];
meanPos = nan(size(elementData(:,1),1),3);%[];
% angXZ = nan(size(elementData(:,1)));%[];
cMap = jet(128);
% Cartesian axes
yAx = [0 1 0]';
% zAx = [0 0 1]';
%

fiberLengthProjected = nan(    length(unique(elementData(:,5) ) ),1);
fiberLength = nan(    length(unique(elementData(:,5) ) ),1);
fiberCurl = nan(    length(unique(elementData(:,5) ) ),1);

counter = 1;
for xLoop = 1:sum(newFiberMarker)-1   % Goes through the data fiber by fiber (useful for measuring fiber length)
   startPos = newFiberIndex(xLoop);
   endPos = newFiberIndex(xLoop+1)-1;

   fiberLengthProjected(xLoop) =  sqrt( (nodalData(elementData(endPos,3),2) - nodalData(elementData(startPos,2),2))^2 + ...
                                  (nodalData(elementData(endPos,3),3) - nodalData(elementData(startPos,2),3))^2 + ...
                                  (nodalData(elementData(endPos,3),4) - nodalData(elementData(startPos,2),4))^2 );

    fiberLength(xLoop) = sum(sqrt((nodalData(elementData(startPos:endPos,3),2) - nodalData(elementData(startPos:endPos,2),2)).^2 + ...
                                  (nodalData(elementData(startPos:endPos,3),3) - nodalData(elementData(startPos:endPos,2),3)).^2 + ...
                                  (nodalData(elementData(startPos:endPos,3),4) - nodalData(elementData(startPos:endPos,2),4)).^2));

   if fiberLengthProjected(xLoop) > fiberLength(xLoop)
    
    figure
    plot(nodalData(elementData(startPos:endPos,3),2) , nodalData(elementData(startPos:endPos,3),3),'-o')
    axis equal
    nodalData(elementData(startPos:endPos,3),1)
    nodalData(elementData(startPos:endPos-1,2),1)
    length(sqrt((nodalData(elementData(startPos+1:endPos,3),2) - nodalData(elementData(startPos:endPos-1,2),2)).^2 + ...
                                  (nodalData(elementData(startPos+1:endPos,3),3) - nodalData(elementData(startPos:endPos-1,2),3)).^2 + ...
                                  (nodalData(elementData(startPos+1:endPos,3),4) - nodalData(elementData(startPos:endPos-1,2),4)).^2))
    disp(stop)
   end

   fiberCurl(xLoop) = fiberLength(xLoop) / fiberLengthProjected(xLoop) - 1.0;

   yLoop = -1;

   while yLoop <=(endPos-startPos)-1 % Goes through each fiber, element by element

       yLoop = yLoop + 1;
       ia = startPos+yLoop;

       % This conditional checks that the start and end node of the element
       % in question is inside the region of interest so that only that
       % material is added. Since many fibers are partially inside, a while
       % loop is used to incrementally "step through" the fiber elements.
       %
       % This implementation is by no means the fastest, but it is easy
       % enough to understand and debug.
       conditionalOne = nodalData(elementData(ia,2),2) > xDim(1) & nodalData(elementData(ia,2),2) < xDim(2) & ...
                        nodalData(elementData(ia,2),3) > yDim(1) & nodalData(elementData(ia,2),3) < yDim(2) ;
       conditionalTwo = nodalData(elementData(ia,3),2) > xDim(1) & nodalData(elementData(ia,3),2) < xDim(2) & ...
                        nodalData(elementData(ia,3),3) > yDim(1) & nodalData(elementData(ia,3),3) < yDim(2) ;
       conditionalThree = 1;%nodalData(elementData(ia,3),4) < (zDim(1) + 0.25*(zDim(2)-zDim(1)));% & nodalData(elementData(ia,3),2) < zDim(2);

       % If the element is entirely inside the region of interest, add it
       % to the calculation of fiber properties.
       if conditionalOne && conditionalTwo && conditionalThree
          % Structural anisotropy module:
          % Goal: Measure the direction of each element to determine the
          % average in-plane orientation
          %
          % Elements are simplified to linear segments.

          tempRand = rand(1);
          if tempRand > 0.5
             sIdx = 3;
             eIdx = 2;
          else
             sIdx = 2;
             eIdx = 3;
          end

          elementVector = [ nodalData(elementData(ia,eIdx),2) - nodalData(elementData(ia,sIdx),2) ;
                            nodalData(elementData(ia,eIdx),3) - nodalData(elementData(ia,sIdx),3) ;
                            nodalData(elementData(ia,eIdx),4) - nodalData(elementData(ia,sIdx),4) ];
          elementVectorNorm = elementVector./norm(elementVector);
          elementVectorRed = elementVector(1:2)./norm(elementVector(1:2)); % In plane vector

          angYX(counter) = acos(dot([1 0],elementVectorRed)); % pi/2-
          meanPos(counter,:) = mean([ nodalData(elementData(ia,eIdx),2) , nodalData(elementData(ia,sIdx),2) ;
                                 nodalData(elementData(ia,eIdx),3) , nodalData(elementData(ia,sIdx),3) ;
                                 nodalData(elementData(ia,eIdx),4) , nodalData(elementData(ia,sIdx),4) ],2);

          % Find quadrant to look in

          if elementVectorRed(1) > 0 && elementVectorRed(2) > 0 % Here we can use any function
              angYX(counter) = acos(elementVectorRed(1));
          elseif elementVectorRed(1) > 0 && elementVectorRed(2) < 0 % Bottom right quadrant
              angYX(counter) = asin(elementVectorRed(2));
          elseif elementVectorRed(1) < 0 && elementVectorRed(2) > 0 % Top left quadrant
              angYX(counter) = acos(elementVectorRed(1));
          elseif elementVectorRed(1) < 0 && elementVectorRed(2) < 0 % Bottom left quadrant, here we need to
              angYX(counter) = pi-asin(elementVectorRed(2));
          elseif elementVectorRed(1) == 0 && elementVectorRed(2) == 1
              angYX(counter) = pi/2;
          elseif elementVectorRed(1) == 1 && elementVectorRed(2) == 0
              angYX(counter) = 0;
          elseif elementVectorRed(1) == 0 && elementVectorRed(2) == -1
              angYX(counter) = -pi/2;
          elseif elementVectorRed(1) == -1 && elementVectorRed(2) == 0
              angYX(counter) = 180;
          else
              disp(stop)
          end

          counter = counter + 1;


            if 0%1 % Debug clause
               subplot(1,3,1)
               plot([nodalData(elementData(ia,2),2) nodalData(elementData(ia,3),2)],[nodalData(elementData(ia,2),3) nodalData(elementData(ia,3),3)],'s--')
               hold on
               plot([-inf inf],[0 0 ],'k--')
               plot([-0 0],[-inf inf],'k--')
               axis equal
               xlim(xDim); ylim(yDim);
               xlabel x; ylabel y;


               subplot(1,3,2)
               plot([0 elementVectorRed(1)],[0 elementVectorRed(2)],'s-r')
               hold on
               plot([0 1],[0 0],'k--','linewidth',2)
               plot(cos(angYX(1:(counter-2))), sin(angYX(1:(counter-2))),'b.')
               title(['Proposed Angle = ' num2str(rad2deg(angYX(counter-1)))])
               axis equal
               xlabel x; ylabel y;
               xlim([-1 1]); ylim([-1 1])
               pause(0.15)
               hold off

               subplot(1,3,3)
               polarhistogram(angYX(1:(counter-1)),360,'displaystyle','bar','Normalization','pdf');


%                if mod(counter-1,100) == 0
%                pause(1)
%                end

            end



        else % If either start nor end is in the zone

       end
   end
end

angYX(isnan(angYX)) = [];
meanPos(isnan(meanPos(:,1)),:) = [];


% Create output structure which is useful for filtering operations
networkStructure.name = networkName;
networkStructure.dir = networkDir;
networkStructure.directionData = angYX;
networkStructure.meanPos = meanPos;

networkStructure.nodalData = nodalData;
networkStructure.realSetData = realSetData;
networkStructure.elementData = elementData;

networkStructure.fiberLengthProjected = fiberLengthProjected;
networkStructure.fiberLength = fiberLength;
networkStructure.fiberCurl = fiberCurl;

% Fiber orientation data
% Note that orientation is only calculated for 1 quadrant and then expanded
% via double-symmetry (A fiber oriented at 0 deg == a fiber oriented at 180
% deg).
% In-plane
networkStructure.angYX = networkStructure.directionData(:,1);
% networkStructure.angYX = [mod(mod(networkStructure.directionData(:,1),2*pi),pi) ; pi+mod(mod(networkStructure.directionData(:,1),2*pi),pi)];% ; ...


