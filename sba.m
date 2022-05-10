function [p, s] = sba(ax, vert_url, faces_url, labels_url, area2highlight, ap_offset, ml_offset, p)
%SBA  Plot Scalable Brain Atlas for macaque
%
% Syntax:
%   p = plot.sba(ax);
%   p = plot.sba(ax, vert_url, faces_url, labels_url);
%   p = plot.sba(__, area2highlight, ap_offset, ml_offset, p);
%
%
% Inputs:
%   ax - (Optional) axes object to plot these things on
%   vert_url - (Optional) char or string -- the location of vertices csv
%   faces_url - (Optional) char or string -- location of faces csv file
%   labels_url - (Optional) char or string -- location of labels csv file
%   area2highlight - (Optional) cell array of char or string array -- One
%                   element for each area to highlight. Default is ["PrG",
%                   "PoG"] (precentral gyrus, postcentral gyrus -- names
%                   depend upon what is used in labels csv file).
%   ap_offset - (Optional; default = 10 mm) -- Offset applied to convert
%               each vertex y-component to ear-bar zero (from original AP
%               "zero").
%   ml_offset - (Optional; default = 0 mm) -- Offset applied to convert
%               each vertex x-component to ear-bar zero (from original ML
%               "zero").
%   p - (Optional) -- If provided, then skips import and just uses data in
%                       patch object p to apply correct formatting.
%
% Output:
%   p - Patch object generated by plotting.
%   s - "Structified" version of patch object.
%
% See also: Contents, example_SBA

if nargin < 1
    ax = gca; 
end

if nargin < 2 || isempty(vert_url)
    vert_url = 'https://scalablebrainatlas.incf.org//templates/DB09/meshes/wholebrain_vertices.csv'; 
end

if nargin < 3 || isempty(faces_url)
    faces_url = 'https://scalablebrainatlas.incf.org//templates/DB09/meshes/wholebrain_faces.csv';
end

if nargin < 4 || isempty(labels_url)
    labels_url = 'https://scalablebrainatlas.incf.org//templates/DB09/meshdata/wholebrain_labels.csv';
end

if nargin < 5 || isempty(area2highlight)
    area2highlight = ["PrG", "PoG"]; 
end

if nargin < 6 || isempty(ap_offset) || isnan(ap_offset)
    ap_offset = 20; 
end

if nargin < 7 || isempty(ml_offset) || isnan(ml_offset)
    ml_offset = 0; 
end

if nargin > 7
    faces = p.Faces;
    vertices = p.Vertices;
    vtx2id = p.FaceVertexCData;
    cmapdata = p.UserData;
else
    % % % Import data from scalable brain atlas (or other) % % %
    options = weboptions('ContentReader', @importdata);
    vertices = webread(vert_url,options);
    faces = webread(faces_url,options);
    vtx2id = webread(labels_url,options);
    vertices(:, 1) = vertices(:, 1) + ml_offset;
    vertices(:, 2) = vertices(:, 2) + ap_offset;
    vertices(:, 3) = vertices(:, 3) - max(vertices(:, 3));
    faces = round(faces) + 1; 
    % % % Set up indexing etc. in MATLAB % % %
    T = getfield(load('Macaque-Brain-Regions.mat', 'T'), 'T');
    id2use = [];
    area2highlight = string(area2highlight);
    for ii = 1:numel(area2highlight)
         id2use = [id2use, T.Index(T.Area == area2highlight(ii))]; %#ok<AGROW>
    end
    
    % % % Set up labeling % % %
    vtx2id(~ismember(vtx2id, id2use)) = 1;
    for ii = 1:numel(id2use)
        vtx2id(vtx2id == id2use(ii)) = (ii + 1);
    end
    cmapdata = []; % We don't have a preset cmapdata from .UserData field
end

% % % Create graphics objects % % %
fig = ax.Parent;
while ~isa(fig, 'matlab.ui.Figure')
    fig = fig.Parent;
end
set(ax, 'NextPlot', 'add');
p = patch(ax, ...
    'Vertices',vertices, ...
    'Faces',faces,...
    'FaceVertexCData',vtx2id,...
    'FaceColor','flat','FaceLighting','phong', ...
    'EdgeColor','none','CDataMapping','scaled');
material(p, 'metal');


xlabel(ax, 'ML (mm)', 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 14);
ylabel(ax, 'AP (mm)', 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 14);
zlabel(ax, 'DV (mm)', 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 14);



if isempty(cmapdata)
    if numel(area2highlight) == 2
        cmapdata = [0.5 0.5 0.5; 
                    1.0 0.0 1.0; 
                    0.0 1.0 1.0];
    else
        cmapdata = [0.5 0.5 0.5; rand(numel(area2highlight), 3)];
    end
end
p.UserData = cmapdata;
colormap(ax, cmapdata);
for ii = 2:size(cmapdata, 1)
    text(ax, -2 - 15*(ii-2), -5-4*(ii-2), 1 - ii, area2highlight(ii-1),...
        'FontName', 'Tahoma', ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'Color', cmapdata(ii, :));
end

lighting(ax, 'gouraud');
view(ax, [-50 23]);
camlight(ax, -50, 23, 'infinite');
xlim(ax, [ax.XLim(1), 5]);
zlim(ax, [ax.ZLim(1), 20]);
daspect(ax, [1 1 1]);

set(fig,'renderer','opengl');
if nargout > 1
    s = utils.sba_patch_2_struct(p); 
end
end