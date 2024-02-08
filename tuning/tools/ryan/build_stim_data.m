function params = build_stim_data(params)

switch size(params.param, 2)
    case 2
        params.p1 = unique(param(:, 1)); %#ok<*USENS>
        params.p2 = unique(param(:, 2));
        params.np1 = numel(params.p1);
        params.np2 = numel(params.p2);
        [params.pp1, params.pp2] = meshgrid(params.p1, params.p2);
    case 3
        params.p1 = unique(param(:, 1)); %#ok<*USENS>
        params.p2 = unique(param(:, 2));
        params.p3 = unique(param(:, 3));
        params.np1 = numel(params.p1);
        params.np2 = numel(params.p2);
        params.np3 = numel(params.p3);
        [params.pp1, params.pp2, params.pp3] = meshgrid(params.p1, params.p2, params.p3);
    case 4
        params.p1 = unique(param(:, 1)); %#ok<*USENS>
        params.p2 = unique(param(:, 2));
        params.p3 = unique(param(:, 3));
        params.p4 = unique(param(:, 4));
        params.np1 = numel(params.p1);
        params.np2 = numel(params.p2);
        params.np3 = numel(params.p3);
        params.np4 = numel(params.p4);
        [params.pp1, params.pp2, params.pp3, params.pp4] = ndgrid(params.p1, params.p2, params.p3, params.p4);
end
end