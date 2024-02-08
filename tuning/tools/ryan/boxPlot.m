
function [xCoordinates,lgdObject] = boxPlot(inputData,NameValueArgs)
% BOXPLOT Plot boxplot from vector or matrix input.
%
% MIT License
% Copyright (c) 2022 Ryan Gorzek
% https://github.com/gorzek-ryan/matlab_viz/blob/main/LICENSE
% https://opensource.org/licenses/MIT
%
% Dependencies: none
% 
% Input Arguments:
%
%   inputData -- vector or matrix of data for boxplot. If inputData is a matrix, 
%                each column will be plotted as a box.If inputData is a vector, specify
%                categorical name-value argument inputLabels to plot multiple boxes. 
%
%   Name-Value Arguments (default):
%
%     inputLabels ([]) -- vector of categorical labels for vector input.
%
%     groupSize (1) -- scalar that specifies the number of boxes by which to group the data (if grouping).
%
%     labelGroup (false) -- logical that specifies whether to shrink the number of x-axis labels to one per group (if grouping).
%
%     boxLabels ([1:nBoxes]) -- cell array of strings that specify x-axis labels.
%
%     boxColors ([0.7, 0.7, 0.7]) -- cell array of RGB vectors that specify box colors.
%
%     boxEdgeColors ([0.0, 0.0, 0.0]) -- 
%
%     boxLineWidth (1) -- scalar that specifies the line width of the box edges and whiskers.
%
%     whiskerColors ([0.0, 0.0, 0.0]) -- 
%
%     outlierSize (30) -- scalar that specifies the size of outlier points.
%
%     plotPoints (false) -- 
%
%     pointSize (10) -- 
%
%     pointColor ([0.0, 0.0, 0.0]) -- 
%
%     connectGroups (false) -- 
%
%     connectLineWidth () -- 
%
%     connectLineColor ([0.0, 0.0, 0.0]) -- 
%
%     plotLegend (false) -- logical that specifies whether or not to plot a legend with the specified parameters.
%
%     lgdLabels (none) -- cell array of strings that specify legend labels.
%
%     lgdColors (none) -- cell array of RBG vectors that specify legend colors.
%
%     lgdColumns (1) -- scalar that specifies the number of columns in the legend.
%
%     lgdOrientation ("horizontal") -- 
%
%     lgdBox ("off") -- string that specifies whether to show legend box. Options are "on" or "off".
%
%     lgdFontSize (12) -- scalar that specifies legend font size.
%
%     lgdLineHeight (1) -- scalar that specifies the height of legend markers.
%
%     lgdLineWidth (1) -- scalar that specifies the width of legend markers.
%
%     lgdLocation ("northeast") -- string that specifies legend location. See MATLAB legend documentation for options.
%
%     lgdPosition ([]) -- 
%
% Output Arguments:
%
%   xCoors -- vector of x-axis coordinates corresponding to each plotted box.
%
%   lgdObject -- legend object. See MATLAB legend documentation for more information.
%

arguments
    
    inputData (:,:) {mustBeNumeric}
    
    NameValueArgs.inputLabels {mustBeVector, mustBeNumeric} = reshape(repmat(1:size(inputData,2),[size(inputData,1),1]),[],1)

    NameValueArgs.groupSize (1,1) {mustBeNumeric} = 1
    NameValueArgs.labelGroups (1,1) logical = false

    NameValueArgs.boxLabels {mustBeA(NameValueArgs.boxLabels,"cell")} = {}
    NameValueArgs.boxColors {mustBeA(NameValueArgs.boxColors,"cell")} = {}
    NameValueArgs.boxAlpha (1,1) {mustBeInRange(NameValueArgs.boxAlpha,0,1)} = 1
    NameValueArgs.boxEdgeColors {mustBeA(NameValueArgs.boxEdgeColors,"cell")} = {}
    NameValueArgs.boxEdgeWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.boxEdgeStyle (1,1) string = "-"
    NameValueArgs.boxEdgeAlpha (1,1) {mustBeInRange(NameValueArgs.boxEdgeAlpha,0,1)} = 1
    NameValueArgs.boxSpacing (1,1) {mustBeNumeric} = 1
    NameValueArgs.boxCurvature (1,2) double = [0,0]

    NameValueArgs.medianColors (1,:) {mustBeA(NameValueArgs.medianColors,"cell")} = {}
    NameValueArgs.medianWidth (1,1) {mustBeNumeric} = 2
    NameValueArgs.medianStyle (1,1) string = "-"
    NameValueArgs.medianAlpha (1,1) {mustBeInRange(NameValueArgs.medianAlpha,0,1)} = 1

    NameValueArgs.whiskerColors (1,:) {mustBeA(NameValueArgs.whiskerColors,"cell")} = {}
    NameValueArgs.whiskerWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.whiskerStyle (1,1) string = "-"
    NameValueArgs.whiskerAlpha (1,1) {mustBeInRange(NameValueArgs.whiskerAlpha,0,1)} = 1

    NameValueArgs.outlierColors (1,:) {mustBeA(NameValueArgs.outlierColors,"cell")} = {}
    NameValueArgs.outlierSize (1,1) {mustBeNumeric} = 30
    NameValueArgs.outlierStyle (1,1) string = "o"
    NameValueArgs.outlierAlpha (1,1) {mustBeInRange(NameValueArgs.outlierAlpha,0,1)} = 1
    NameValueArgs.outlierJitter {mustBeMember(NameValueArgs.outlierJitter,["none","rand","randn"])} = "none"

    NameValueArgs.plotPoints {mustBeNumericOrLogical} = false
    NameValueArgs.pointColors {mustBeA(NameValueArgs.pointColors,"cell")} = {}
    NameValueArgs.pointSize (1,1) {mustBeNumeric} = 20
    NameValueArgs.pointStyle (1,1) string = "."
    NameValueArgs.pointAlpha (1,1) {mustBeInRange(NameValueArgs.pointAlpha,0,1)} = 1
    NameValueArgs.pointJitter {mustBeMember(NameValueArgs.pointJitter,["none","rand","randn"])} = "none"

    NameValueArgs.jitterBound (1,1) {mustBeNumeric} = 0.75

    NameValueArgs.plotLines {mustBeNumericOrLogical} = false
    NameValueArgs.lineColors {mustBeA(NameValueArgs.lineColors,"cell")} = {}
    NameValueArgs.lineWidth (1,1) {mustBeNumeric} = 0.5
    NameValueArgs.lineStyle (1,1) string = "-"
    NameValueArgs.lineAlpha (1,1) {mustBeInRange(NameValueArgs.lineAlpha,0,1)} = 1

    NameValueArgs.plotLegend (1,1) logical = false
    NameValueArgs.lgdLabels (1,:) {mustBeA(NameValueArgs.lgdLabels,"cell")} = {}
    NameValueArgs.lgdColors (1,:) {mustBeA(NameValueArgs.lgdColors,"cell")} = {}
    NameValueArgs.lgdColumns (1,1) {mustBeNumeric} = 1
    NameValueArgs.lgdOrientation {mustBeMember(NameValueArgs.lgdOrientation,["vertical","horizontal"])} = "horizontal"
    NameValueArgs.lgdBox {mustBeMember(NameValueArgs.lgdBox,["on","off"])} = "off"
    NameValueArgs.lgdFontSize (1,1) {mustBeNumeric} = 12
    NameValueArgs.lgdLineHeight (1,1) {mustBeNumeric} = 1
    NameValueArgs.lgdLineWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.lgdLocation (1,1) string = "northeast"
    NameValueArgs.lgdPosition (1,4) {mustBeInRange(NameValueArgs.lgdPosition,0,1)} = [0,0,0,0]

end

% Get number of boxes, groups, samples (including NaN), and missing values.
nBoxes = numel(unique(NameValueArgs.inputLabels));
nGroups = nBoxes/NameValueArgs.groupSize;
nSamples = size(inputData,1);
nMissing = nnz(all(isnan(inputData),2));

% Throw error if number of groups is not an integer.
if rem(nGroups,1) ~= 0
    error("Number of input boxes is not divisible by number of groups."); 
end

% Reshape input data into a column vector and get unique labels.
inputData = reshape(inputData,[],1);
uniqueLabels = sort(unique(NameValueArgs.inputLabels),"ascend")';

% Set default box labels (numbered) if none are specified.
if isempty(NameValueArgs.boxLabels) && ...
   NameValueArgs.groupSize == 1

    NameValueArgs.boxLabels = cellstr(string(1:nBoxes));

elseif isempty(NameValueArgs.boxLabels) && ...
       NameValueArgs.groupSize > 1 && ...
       NameValueArgs.labelGroups == false

    NameValueArgs.boxLabels = cellstr(string(1:NameValueArgs.groupSize));

elseif isempty(NameValueArgs.boxLabels) && ...
       NameValueArgs.labelGroups == true
    
    NameValueArgs.boxLabels = cellstr(string(1:nGroups));

end

% Set default box colors (gray) if none are specified, or replicate
% single or group-level color specification.
if isempty(NameValueArgs.boxColors)
    NameValueArgs.boxColors = repmat({[0.7,0.7,0.7]},[1,nBoxes]);
elseif numel(NameValueArgs.boxColors) == 1
    NameValueArgs.boxColors = repmat(NameValueArgs.boxColors,[1,nBoxes]);
elseif numel(NameValueArgs.boxColors) == NameValueArgs.groupSize
    NameValueArgs.boxColors = repmat(NameValueArgs.boxColors,[1,nGroups]);
end

% Set default box edge colors (black) if none are specified, or replicate
% single or group-level color specification.
if isempty(NameValueArgs.boxEdgeColors)
    NameValueArgs.boxEdgeColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NameValueArgs.boxEdgeColors) == 1
    NameValueArgs.boxEdgeColors = repmat(NameValueArgs.boxEdgeColors,[1,nBoxes]);
elseif numel(NameValueArgs.boxEdgeColors) == NameValueArgs.groupSize
    NameValueArgs.boxEdgeColors = repmat(NameValueArgs.boxEdgeColors,[1,nGroups]);
end

% Set default whisker colors (black) if none are specified, or replicate
% single or group-level color specification.
if isempty(NameValueArgs.whiskerColors)
    NameValueArgs.whiskerColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NameValueArgs.whiskerColors) == 1
    NameValueArgs.whiskerColors = repmat(NameValueArgs.whiskerColors,[1,nBoxes]);
elseif numel(NameValueArgs.whiskerColors) == NameValueArgs.groupSize
    NameValueArgs.whiskerColors = repmat(NameValueArgs.whiskerColors,[1,nGroups]);
end

% Set default median colors (black) if none are specified, or replicate
% single or group-level color specification.
if isempty(NameValueArgs.medianColors)
    NameValueArgs.medianColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NameValueArgs.medianColors) == 1
    NameValueArgs.medianColors = repmat(NameValueArgs.medianColors,[1,nBoxes]);
elseif numel(NameValueArgs.medianColors) == NameValueArgs.groupSize
    NameValueArgs.medianColors = repmat(NameValueArgs.medianColors,[1,nGroups]);
end

% Set default outlier colors (matched to box colors) if none are specified,
% or replicate single or group-level color specification.
if isempty(NameValueArgs.outlierColors)
    NameValueArgs.outlierColors = NameValueArgs.boxColors;
elseif numel(NameValueArgs.outlierColors) == 1
    NameValueArgs.outlierColors = repmat(NameValueArgs.outlierColors,[1,nBoxes]);
elseif numel(NameValueArgs.outlierColors) == NameValueArgs.groupSize
    NameValueArgs.outlierColors = repmat(NameValueArgs.outlierColors,[1,nGroups]);
end

% Set default point colors (black) if not specified.
if isempty(NameValueArgs.pointColors)
    NameValueArgs.pointColors = repmat({[0.0,0.0,0.0]},[nSamples,nBoxes]);
elseif numel(NameValueArgs.pointColors) == 1
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors,[nSamples,nBoxes]);
elseif isvector(NameValueArgs.pointColors) && ...
       any(size(NameValueArgs.pointColors) == NameValueArgs.groupSize)
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors,[nSamples,nGroups]);
elseif isvector(NameValueArgs.pointColors) && ...
       any(size(NameValueArgs.pointColors) == nBoxes)
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors,[nSamples,1]);
elseif all(ismember(size(NameValueArgs.pointColors),[nSamples - nMissing,NameValueArgs.groupSize]))
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors,[1,nGroups]);
end

% Check default color stream specification for point colors.
point_defaultIdx = cell2mat(cellfun(@(x) strcmp(x,"default"),NameValueArgs.pointColors,"UniformOutput",false));
if any(point_defaultIdx,"all")
    defaultColors = repmat(num2cell(colororder,2),[ceil(size(NameValueArgs.pointColors,1)/7), ...
                                                   size(NameValueArgs.pointColors,2)]);
    defaultColors = defaultColors(1:size(NameValueArgs.pointColors,1),:);
    NameValueArgs.pointColors(point_defaultIdx) = defaultColors(point_defaultIdx);
end

% Set default line colors (black) if not specified.
nLines = nGroups*(NameValueArgs.groupSize - 1);
if isempty(NameValueArgs.lineColors)
    NameValueArgs.lineColors = repmat({[0.0,0.0,0.0]},[nSamples,nLines]);
elseif numel(NameValueArgs.lineColors) == 1
    NameValueArgs.lineColors = repmat(NameValueArgs.lineColors,[nSamples,nLines]);
elseif isvector(NameValueArgs.lineColors) && ...
       any(size(NameValueArgs.lineColors) == nGroups)
    NameValueArgs.lineColors = NameValueArgs.lineColors(kron(1:nGroups,ones(nSamples,groupSize-1)));
elseif isvector(NameValueArgs.lineColors) && ...
       any(size(NameValueArgs.lineColors) == nLines)
    NameValueArgs.lineColors = repmat(NameValueArgs.lineColors,[nSamples,1]);
elseif all(ismember(size(NameValueArgs.lineColors),[nSamples - nMissing,nGroups]))
    NameValueArgs.lineColors = NameValueArgs.lineColors(:,kron(1:nGroups,ones(1,groupSize-1)));
end

% Check default color stream specification for line colors.
line_defaultIdx = cell2mat(cellfun(@(x) strcmp(x,"default"),NameValueArgs.lineColors,"UniformOutput",false));
if any(line_defaultIdx,"all")
    defaultColors = repmat(num2cell(colororder,2),[ceil(size(NameValueArgs.lineColors,1)/7), ...
                                                   size(NameValueArgs.lineColors,2)]);
    defaultColors = defaultColors(1:size(NameValueArgs.lineColors,1),:);
    NameValueArgs.lineColors(line_defaultIdx) = defaultColors(line_defaultIdx);
end

% Check boxLabels against nGroups to warn user about specifying labelGroups.
if NameValueArgs.labelGroups == true && ...
   numel(NameValueArgs.boxLabels) > nGroups
    error(["Number of box labels exceeds number of groups, " ...
           "did you mean to specify labelGroups = false?"]);
elseif NameValueArgs.labelGroups == true && ...
       numel(NameValueArgs.boxLabels) < nGroups
    error("Insufficient number of box labels for number of groups.");
elseif NameValueArgs.labelGroups == false && ...
       NameValueArgs.groupSize ~= 1 && ...
       numel(NameValueArgs.boxLabels) == nGroups && ...
       numel(NameValueArgs.boxColors) == NameValueArgs.groupSize
    error(["Insufficient number of box labels for number of groups, " ...
           "did you mean to specify labelGroups = true?"]);
end

% Generate x-coordinates for plotting boxes.
if nGroups == 1
    xCoordinates = (0.70 : 0.65 : 0.70 + (0.65*nBoxes) - 0.65).*boxSpacing;
else
    xCoordinates = 0.70 : 0.55 : 0.70 + (0.55*(nBoxes/nGroups) - 0.55); 
    initCoors = 0.70 : 0.55 : 0.70 + (0.55*(nBoxes/nGroups) - 0.55);
    for grp = 2:nGroups
        xCoordinates = horzcat(xCoordinates, initCoors + xCoordinates(end) + 0.4);
    end
    xCoordinates = xCoordinates.*NameValueArgs.boxSpacing;
end

% Generate matrix of jitter if specified.
if strcmp(NameValueArgs.outlierJitter,"rand") || ...
   strcmp(NameValueArgs.pointJitter,"rand")

    jitterMat = rand(size(reshape(inputData,[],nBoxes))).*(0.5*NameValueArgs.jitterBound);
    jitterMat = jitterMat - mean(jitterMat,1);

elseif strcmp(NameValueArgs.outlierJitter,"randn") || ...
       strcmp(NameValueArgs.pointJitter,"randn")

    randnMat = randn(size(reshape(inputData,[],nBoxes)));
    jitterMat = (randnMat./max(randnMat,[],1)).*(0.5*NameValueArgs.jitterBound);
    jitterMat = jitterMat - mean(jitterMat,1);

elseif strcmp(NameValueArgs.outlierJitter,"none") && ...
       strcmp(NameValueArgs.pointJitter,"none")

    jitterMat = zeros(size(reshape(inputData,[],nBoxes)));

end

% Intialize matrices for storing max/min values to set axis limits.
maxMat = zeros(nBoxes,3); minMat = zeros(nBoxes,3);

% Plot boxes.
for box = uniqueLabels

    % Get current box location and data.
    boxNum = find(uniqueLabels == box);
    currData = inputData(NameValueArgs.inputLabels == box,1);
    
    % Plot box if there are at least 4 data points.
    if nnz(~isnan(currData)) > 4

        boxMedian = median(currData,1,"omitnan");

        upperQuantile = quantile(currData, 0.75);
        lowerQuantile = quantile(currData, 0.25);

        maxWhisker = upperQuantile + 1.5*(upperQuantile - lowerQuantile); 
        upperWhisker = max(currData(currData < maxWhisker & currData >= upperQuantile));
        if isempty(upperWhisker), upperWhisker = maxWhisker; end

        minWhisker = lowerQuantile - 1.5*(upperQuantile - lowerQuantile);
        lowerWhisker = min(currData(currData > minWhisker & currData <= lowerQuantile));
        if isempty(lowerWhisker), lowerWhisker = minWhisker; end
        
        hold on;

        % Plot box with bounds at quartiles.
        rectangle("Position",  [xCoordinates(boxNum)-0.25, lowerQuantile, 0.5, upperQuantile-lowerQuantile],...
                  "FaceColor", [NameValueArgs.boxColors{boxNum}, NameValueArgs.boxAlpha],...
                  "EdgeColor", [NameValueArgs.boxEdgeColors{boxNum}, NameValueArgs.boxEdgeAlpha],...
                  "LineWidth",  NameValueArgs.boxEdgeWidth,...
                  "LineStyle",  NameValueArgs.boxEdgeStyle,...
                  "Curvature",  NameValueArgs.boxCurvature,...
                  "Tag",       "Box");

        % Plot median line.
        line([xCoordinates(boxNum)-0.25, xCoordinates(boxNum)+0.25],...
             [boxMedian, boxMedian],...
             "Color",    [NameValueArgs.medianColors{boxNum}, NameValueArgs.medianAlpha],...
             "LineWidth", NameValueArgs.medianWidth,...
             "LineStyle", NameValueArgs.medianStyle,...
             "Tag",      "Median");

        % Plot lower whisker.
        line([xCoordinates(boxNum), xCoordinates(boxNum)],...
             [lowerQuantile, lowerWhisker],...
             "Color",    [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha],...
             "LineWidth", NameValueArgs.whiskerWidth,...
             "LineStyle", NameValueArgs.whiskerStyle,...
             "Tag",      "Lower Whisker");

        % Plot upper whisker.
        line([xCoordinates(boxNum), xCoordinates(boxNum)],...
             [upperQuantile, upperWhisker],...
             "Color",    [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha],...
             "LineWidth", NameValueArgs.whiskerWidth,...
             "LineStyle", NameValueArgs.whiskerStyle,...
             "Tag",      "Upper Whisker");

        % Plot lower whisker bar.
        line([xCoordinates(boxNum)-0.1, xCoordinates(boxNum)+0.1],...
             [lowerWhisker, lowerWhisker],...
             "Color",    [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha],...
             "LineWidth", NameValueArgs.whiskerWidth,...
             "LineStyle", NameValueArgs.whiskerStyle,...
             "Tag",      "Lower Whisker Bar");

        % Plot upper whisker bar.
        line([xCoordinates(boxNum)-0.1, xCoordinates(boxNum)+0.1],...
             [upperWhisker, upperWhisker],...
             "Color",    [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha],...
             "LineWidth", NameValueArgs.whiskerWidth,...
             "LineStyle", NameValueArgs.whiskerStyle,...
             "Tag",      "Upper Whisker Bar");

        % Plot upper outliers.
        if any(currData > upperWhisker)
            scatter(xCoordinates(boxNum) - jitterMat(currData > upperWhisker,boxNum),...
                    currData(currData > upperWhisker),...
                    NameValueArgs.outlierSize,...
                    "MarkerFaceColor", NameValueArgs.outlierColors{boxNum},...
                    "MarkerEdgeColor", NameValueArgs.outlierColors{boxNum},...
                    "Marker",          NameValueArgs.outlierStyle,...
                    "MarkerFaceAlpha", NameValueArgs.outlierAlpha,...
                    "MarkerEdgeAlpha", NameValueArgs.outlierAlpha,...
                    "Tag",            "Outlier");
            upperOutliers = currData(currData > upperWhisker);
        else
            upperOutliers = nan;
        end
        
        % Plot lower outliers.
        if any(currData < lowerWhisker)

            scatter(xCoordinates(boxNum) - jitterMat(currData < lowerWhisker,boxNum),...
                    currData(currData < lowerWhisker),...
                    NameValueArgs.outlierSize,...
                    "MarkerFaceColor", NameValueArgs.outlierColors{boxNum},...
                    "MarkerEdgeColor", NameValueArgs.outlierColors{boxNum},...
                    "Marker",          NameValueArgs.outlierStyle,...
                    "MarkerFaceAlpha", NameValueArgs.outlierAlpha,...
                    "MarkerEdgeAlpha", NameValueArgs.outlierAlpha,...
                    "Tag",            "Outlier");

            lowerOutliers = currData(currData < lowerWhisker);
        else
            lowerOutliers = nan;
        end

        % Store min/max from each box to set axis limits.
        maxMat(boxNum,:) = [upperQuantile, upperWhisker, max(upperOutliers)];
        minMat(boxNum,:) = [lowerQuantile, lowerWhisker, min(lowerOutliers)];
        
    elseif nnz(~isnan(currData)) == 0
        
        % Store min/max from each box to set axis limits.
        maxMat(boxNum,:) = [repmat(max(currData),[1,3])];
        minMat(boxNum,:) = [repmat(min(currData),[1,3])];
    
    elseif nnz(~isnan(currData)) <= 4

        boxMedian = median(currData,1,"omitnan");
        
        hold on;

        % Plot median line.
        line([xCoordinates(boxNum)-0.25, xCoordinates(boxNum)+0.25],...
             [boxMedian, boxMedian],...
             "Color",    [NameValueArgs.medianColors{boxNum}, NameValueArgs.medianAlpha],...
             "LineWidth", NameValueArgs.medianWidth,...
             "LineStyle", NameValueArgs.medianStyle,...
             "Tag",      "Median");

        % Scatter outliers.
        scatter(xCoordinates(boxNum) - jitterMat(:,boxNum),...
                currData,...
                outlierSize,...
                "MarkerFaceColor", NameValueArgs.outlierColors{boxNum},...
                "MarkerEdgeColor", NameValueArgs.outlierColors{boxNum},...
                "Marker",          NameValueArgs.outlierStyle,...
                "MarkerFaceAlpha", NameValueArgs.outlierAlpha,...
                "MarkerEdgeAlpha", NameValueArgs.outlierAlpha,...
                "Tag",            "Outlier");

        % Store min/max from each box to set axis limits.
        maxMat(boxNum,:) = [repmat(max(currData), [1,3])];
        minMat(boxNum,:) = [repmat(min(currData), [1,3])];
        
    end
    
end

% Connect groups with lines if specified.
if (isnumeric(NameValueArgs.plotLines) || NameValueArgs.plotLines == true) && ...
   nGroups == nBoxes

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(NameValueArgs.plotLines)
        lineIdx = NameValueArgs.plotLines; 
    else
        lineIdx = 1:nBoxes; 
    end

    for connection = 1:numel(lineIdx)-1
        for sample = find(~all(isnan(currData),2))'

            plot(xCoordinates(lineIdx(connection:connection+1)) - jitterMat(sample,lineIdx(connection:connection+1)),...
                 currData(sample,lineIdx(connection:connection+1)),...
                 "Color",    [NameValueArgs.lineColors{find(~all(isnan(currData),2)) == sample,connection}, NameValueArgs.lineAlpha],...
                 "LineWidth", NameValueArgs.lineWidth,...
                 "LineStyle", NameValueArgs.lineStyle,...
                 "Tag",      "Line");

        end
    end

elseif (isnumeric(NameValueArgs.plotLines) || NameValueArgs.plotLines == true) && ...
       nGroups ~= nBoxes

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(NameValueArgs.plotLines)
        lineIdx = NameValueArgs.plotLines;
    else
        lineIdx = 1:nBoxes/nGroups;
    end

    connectionNum = 1;
    for group = 1:nGroups

        boxIdx = reshape(1:nBoxes,[],nGroups)'; 
        boxIdx = boxIdx(:,lineIdx);

        for connection = 1:numel(boxIdx(group,:))-1
            for sample = find(~all(isnan(currData),2))'

                plot(xCoordinates(boxIdx(group,connection:connection + 1)) - jitterMat(sample,boxIdx(group,connection:connection + 1)),...
                     currData(sample,boxIdx(group,connection:connection + 1)),...
                     "Color",    [NameValueArgs.lineColors{find(~all(isnan(currData),2)) == sample,connectionNum}, NameValueArgs.lineAlpha],...
                     "LineWidth", NameValueArgs.lineWidth,...
                     "LineStyle", NameValueArgs.lineStyle,...
                     "Tag",      "Line");
            
            end
            connectionNum = connectionNum + 1;
        end

    end

end

% Plot individual data points if specified.
if (isnumeric(NameValueArgs.plotPoints) || NameValueArgs.plotPoints == true)

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(NameValueArgs.plotPoints)
        pointIdx = NameValueArgs.plotPoints;
    else
        pointIdx = 1:nBoxes;
    end

    for sample = find(~all(isnan(currData),2))'

        scatter(xCoordinates(pointIdx) - jitterMat(sample,pointIdx),...
                currData(sample,pointIdx),...
                NameValueArgs.pointSize,...
                vertcat(NameValueArgs.pointColors{find(~all(isnan(currData),2)) == sample,pointIdx}),...
                "Marker",          NameValueArgs.pointStyle,...
                "MarkerFaceAlpha", NameValueArgs.pointAlpha,...
                "MarkerEdgeAlpha", NameValueArgs.pointAlpha,...
                "Tag",            "Point");

    end

end

% Set x-axis limits.
xlim([0.2*NameValueArgs.boxSpacing, xCoordinates(end)+(0.5*NameValueArgs.boxSpacing)]);

% Set x-axis tick labels.
xLabPos = [];
if NameValueArgs.labelGroups == true
    for group = 1:(nBoxes/nGroups):nBoxes
        xLabPos = horzcat(xLabPos, median(xCoordinates(group : group+(nBoxes/nGroups)-1)));
    end
    set(gca,"xtick",      xLabPos,...
            "xticklabel", NameValueArgs.boxLabels);
else 
    set(gca,"xtick",      xCoordinates,...
            "xticklabel", NameValueArgs.boxLabels);
end

% Set figure and axes appearance.
set(gcf,"color",      "w");
set(gca,"box",        "off",...
        "XColor",     "k",...
        "YColor",     "k",...
        "TickDir",    "out",...
        "TickLength", [0.01,0.01],...
        "FontSize",    13,...
        "LineWidth",   1);

% Set y-axis limits.
if ~all(isnan(maxMat),"all")
    yUpper = max(maxMat,[],"all"); yLower = min(minMat,[],"all");
    yExt = (yUpper - yLower)*0.2;
    ylim([yLower-yExt, yUpper+yExt]);
end

% Plot legend if specified.
if NameValueArgs.plotLegend == true && ...
   ~isempty(NameValueArgs.lgdLabels) && ...
   ~isempty(NameValueArgs.lgdColors)

    if NameValueArgs.lgdFontSize > 2
        NameValueArgs.lgdLineHeight = (NameValueArgs.lgdFontSize - 2) * NameValueArgs.lgdLineHeight;
    end

    for lgdEntry = 1:numel(NameValueArgs.lgdColors)

        currData = inputData(NameValueArgs.inputLabels == lgdEntry,1);
        med = median(currData,1,"omitnan");

        % Plot line with box color for legend.
        line([xCoordinates(lgdEntry)-0.25, xCoordinates(lgdEntry)+0.25],...
             [med, med],...
             "Color",     NameValueArgs.lgdColors{lgdEntry},...
             "LineWidth", NameValueArgs.lgdLineHeight,...
             "Tag",      "Legend Line");
    
        % Plot white line over line with box color for legend.
        line([xCoordinates(lgdEntry)-0.25, xCoordinates(lgdEntry)+0.25],...
             [med, med],...
             "Color",    [1.0,1.0,1.0],...
             "LineWidth", NameValueArgs.lgdLineHeight,...
             "Tag",      "Legend Line Cover");

        set(gca, "Children",circshift(gca().Children,-2,1));

    end

    warning("off","MATLAB:handle_graphics:exceptions:SceneNode");

    lgdObject = legend(findobj(gca, "Tag","Legend Line"),...
                       strcat('\fontsize{',num2str(NameValueArgs.lgdFontSize),'}',NameValueArgs.lgdLabels),...
                       "AutoUpdate", "off",...
                       "NumColumns",  NameValueArgs.lgdColumns,...
                       "Orientation", NameValueArgs.lgdOrientation,...
                       "Box",         NameValueArgs.lgdBox,...
                       "Location",    NameValueArgs.lgdLocation);

    lgdObject.ItemTokenSize = [30 * NameValueArgs.lgdLineWidth, 9];

    if any(NameValueArgs.lgdPosition)
        lgdObject.Position = NameValueArgs.lgdPosition;
    end

else
    
    lgdObject = [];

end

end
