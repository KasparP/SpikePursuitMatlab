function rois = loadRois(filename)
% Load ROIs from .mat file 
%filename = 'Y:\users\Amrita\STVoltron\403106\FOV6_metadata.mat';

s = load(filename); % Structure with rois 
fields = fieldnames(s);
indexes = ~cellfun(@isempty, strfind(fieldnames(s), 'Cell'));
cells = fields(indexes);

for cell = length(cells):-1:1
    rois(:, :, cell) = s.(cells{cell});
end
end