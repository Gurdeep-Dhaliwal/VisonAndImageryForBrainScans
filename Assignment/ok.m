...
% Initialize accuracy_table with an additional row for the Region Growing method
accuracy_table = zeros(5, 10);

...
% Loop through each slice, and apply the segmentation methods
for i = 1:10
    ...
    % Preprocessing for the Region Growing method
    img4 = slice;
    img4 = imgaussfilt(img4, 2);
    img4 = imsharpen(img4, 'Amount', 0.5, 'Radius', 1.5);
    
    % Region Growing segmentation
    X = reshape(img4, [], 1);
    [~, C] = kmeans(X, 5);
    seeds = [];
    for c = 1:numel(C)
        [~, index] = min(abs(img4(:) - C(c)));
        [x, y] = ind2sub(size(img4), index);
        seeds = [seeds; x, y];
    end
    region_growing_labels = regionGrowing(img4, seeds, 1000);
    
    % Remove background
    region_growing_labels(slice < bg_threshold) = 0;

    % Region Growing segmentation accuracy
    region_growing_acc = dice(region_growing_labels, double(label_slice));
    accuracy_table(5, i) = mean(region_growing_acc);
end

