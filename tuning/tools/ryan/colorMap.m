
function [cmap, cmap_path] = colorMap(RGBs, NameValueArgs)
% COLORMAP  Custom colormap creation.
%
% MIT License
% Copyright (c) 2022 Ryan Gorzek
% https://github.com/gorzek-ryan/matlab_viz/blob/main/LICENSE
% https://opensource.org/licenses/MIT
%
arguments
    RGBs
    NameValueArgs.cmap_length (1,1) double = 100
    NameValueArgs.save_cmap {mustBeNumericOrLogical} = false
    NameValueArgs.cmap_name {mustBeText} = "mviz_cmap_%i"
    NameValueArgs.cmap_path {mustBeText} = strcat("mviz_cmaps", filesep)
end

% Create colormap.
if iscell(RGBs)
    RGBs = vertcat(RGBs{:});
end
assert(size(RGBs, 2) == 3, ...
           "Input must be m-by-3 matrix or 1-by-n cell array of RGB vectors.");
num_RGBs = size(RGBs, 1);
query_pts = linspace(1, num_RGBs, NameValueArgs.cmap_length);
R_interp = interp1(1:num_RGBs, RGBs(:,1), query_pts, "linear");
G_interp = interp1(1:num_RGBs, RGBs(:,2), query_pts, "linear");
B_interp = interp1(1:num_RGBs, RGBs(:,3), query_pts, "linear");
cmap = [R_interp', G_interp', B_interp'];

% Figure out where to save colormap M-file, if specified.
if NameValueArgs.save_cmap == true
    % Check the path argument, otherwise put it with this function.
    [filepath, name, ext] = fileparts(NameValueArgs.cmap_path);
    assert((name == "" || string(name) == string(NameValueArgs.cmap_name)),  ...
               "Colormap name must match the filename in the colormap path.");
    assert((ext == "" || string(ext) == ".m"),  ...
               "Colormap must be saved as an M-file (.m).");
    this_function_path = fileparts(mfilename('fullpath'));
    if isfolder(filepath)
        % If a valid full or relative (to current folder) is specified, put
        % colormap M-file there.
        full_path = dir(filepath).folder;
    elseif isfolder(fullfile(this_function_path, filepath))
        % If a valid path relative to the location of this function is specified, 
        % put colormap M-file there.
        full_path = fullfile(this_function_path, filepath);
    else
        % Otherwise, put the M-file in a folder called mviz_cmaps below
        % this function.
        full_path = fullfile(this_function_path, strcat("mviz_cmaps", filesep));
        mkdir(full_path);
    end
    addpath(full_path);
    % Avoid overwriting previous colormaps named with default.
    cmap_name = NameValueArgs.cmap_name;
    if contains(cmap_name, "mviz_cmap_")
        dir_names = string({dir(full_path).name});
        default_name_idx = contains(dir_names, "mviz_cmap_");
        default_nums = 1:1000;
        if any(default_name_idx)
            default_names = dir_names(default_name_idx);
            taken_digits = str2double(extract(default_names, digitsPattern));
            default_nums = setdiff(default_nums, taken_digits);
        end
        cmap_name = sprintf(cmap_name, default_nums(1));
    end
    % Create M-file for colormap to be callable by MATLAB's colormap function.
    cmap_path = fullfile(full_path, cmap_name + ".m");
    fileID = fopen(cmap_path, "w"); % TODO
    cmap_as_str = string(sprintf(insertAfter(mat2str(cmap, 3), ";", " ... \n")));
    fprintf(fileID, get_cmap_function_text(cmap_name, cmap_as_str)); % TODO
    fclose(fileID); % TODO
    savepath(cmap_path);
else
    cmap_path = "";
end

    function cmap_function_text = get_cmap_function_text(cmap_name, cmap_as_str)

        template_text = ...
                        "function map = %s(m)" + "\n" + ...
                        "if nargin < 1" + "\n" + ...
                        "\t" + "f = get(groot,'CurrentFigure');" + "\n" + ...
                        "\t" + "if isempty(f)" + "\n" + ...
                        "\t" + "m = size(get(groot,'DefaultFigureColormap'),1);" + "\n" + ...
                        "\t" + "else" + "\n" + ...
                        "\t" + "m = size(f.Colormap,1);" + "\n" + ...
                        "\t" + "end" + "\n" + ...
                        "end" + "\n" + ...
                        "\n" + ...
                        "map = ..." + "\n" + ...
                        "%s" + "\n" + ...
                        "\n" + ...
                        "end" + "\n";
        cmap_function_text = sprintf(template_text, cmap_name, cmap_as_str);
    end

end
