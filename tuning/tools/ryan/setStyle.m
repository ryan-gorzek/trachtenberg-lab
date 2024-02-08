
function setStyle(NameValueArgs)
% setStyle Set axes/figure style, labels, limits, and position.
%
% Author: Ryan Gorzek
%
% Dependencies: none
%
% Input Arguments:
%
%   Name-Value Arguments (default):
%
%     title (none) -- string that specifies plot title.
%
%     xlabel (none) -- string that specifies x-axis label.
%
%     xticks (automatic) -- vector that specifies x-axis tick marks.
%
%     xticklabels (automatic) -- cell array of strings that specify x-axis tick labels.
%
%     xtickangle (automatic) -- scalar that specifies x-axis tick label angle.
%
%     xlim (automatic) -- vector that specifies x-axis limits.
%
%     ylabel (none) -- string that specifies y-axis label.
%
%     yticks (automatic) -- vector that specifies y-axis tick marks.
%
%     yticklabels (automatic) -- cell array of strings that specify y-axis tick labels.
%
%     ytickangle (automatic) -- scalar that specifies y-axis tick label angle.
%
%     ylim (automatic) -- vector that specifies y-axis limits.
%
%     fontSize (13) -- scalar that specifies font size of axes text.
%
%     lineWidth (1) -- scalar that specifies axes line width.
%
%     figPosition (automatic) -- vector that specifies figure position.
%

arguments
    
    NameValueArgs.title (1,1) string = ''
    NameValueArgs.xlabel (1,1) string = ''
    NameValueArgs.xticks (1,:) double = []
    NameValueArgs.xticklabels = []
    NameValueArgs.xtickangle = []
    NameValueArgs.xlim = []
    NameValueArgs.ylabel (1,1) string = ''
    NameValueArgs.yticks (1,:) double = []
    NameValueArgs.yticklabels = []
    NameValueArgs.ytickangle = []
    NameValueArgs.ylim = []
    NameValueArgs.fontSize (1,1) double = 13
    NameValueArgs.tickLength = [0.01, 0.01]
    NameValueArgs.lineWidth (1,1) double = 1
    NameValueArgs.interpreter (1,1) string = "none"
    NameValueArgs.figPosition = []

end

NVAs = NameValueArgs;

% if ~isempty(NVAs.title), title(NVAs.title, "Interpreter",NVAs.interpreter); end
% if ~isempty(NVAs.xlabel), xlabel(NVAs.xlabel, "Interpreter",NVAs.interpreter); end
% if isnan(NVAs.xticks), xticks([]); elseif ~isempty(NVAs.xticks), xticks(NVAs.xticks); end
% if ~isempty(NVAs.xticklabels), xticklabels(NVAs.xticklabels); end
% if ~isempty(NVAs.xtickangle), xtickangle(NVAs.xtickangle); end
% if ~isempty(NVAs.xlim), xlim(NVAs.xlim); end
% if ~isempty(NVAs.ylabel), ylabel(NVAs.ylabel, "Interpreter",NVAs.interpreter); end
% if isnan(NVAs.yticks), yticks([]); elseif ~isempty(NVAs.yticks), yticks(NVAs.yticks); end
% if ~isempty(NVAs.yticklabels), yticklabels(NVAs.yticklabels); end
% if ~isempty(NVAs.ytickangle), ytickangle(NVAs.ytickangle); end
% if ~isempty(NVAs.ylim), ylim(NVAs.ylim); end

set(gcf,'color','w');
set(gca,'box','off','XColor','k','YColor','k','TickDir','out','TickLength',NVAs.tickLength,'FontSize',NVAs.fontSize,'LineWidth',NVAs.lineWidth);

if ~isempty(NVAs.figPosition), set(gcf,'Position',NVAs.figPosition); end

end
