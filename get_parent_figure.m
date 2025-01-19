function fig = get_parent_figure(h,options)
%GET_PARENT_FIGURE Returns the parent figure even if it's not directly 'Parent' in the tree of the provided graphics object.
%
% Syntax:
%   fig = plot.get_parent_figure(h,'Name',value,...);
%
% Inputs:
%   h - Graphics object (e.g. axes handle, or a line in a plot, etc.)
%   
% Options:
%   'MaxNestingLevel' {mustBeInteger, mustBePositive} = 10 | Increase or decrease to change multiple allowed nesting loop check iterations.
%
% Output:
%   fig - "Parent" figure handle, even if it is not directly 'Parent' of
%           this graphics object (i.e. if axes is in a tiledLayout, etc.)
%
% See also: Contents

arguments
    h (1,1)
    options.MaxNestingLevel {mustBeInteger, mustBePositive} = 10; % Prevent infinite loop
end

fig = h.Parent;
counter = 0;
while ~isa(fig,'matlab.ui.Figure') && (counter < options.MaxNestingLevel)
    fig = fig.Parent;
    counter = counter + 1;
end
if counter == options.MaxNestingLevel
    error("Could not find parent figure at maximum allowed nesting level (%d).", options.MaxNestingLevel);
end

end