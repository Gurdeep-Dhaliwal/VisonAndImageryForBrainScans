load("Brain.mat")


num_clusters = 5;
accuracy_table = zeros(4, 10);
percentages = [0, 7, 19, 36, 58.4, 74.1];

for i = 1:10
    img = T1(:,:,i);
    label_slice = label(:,:,i);
    %% 2D Region Growing
% Normalize image intensities to the range [0, 1]
normalized_img = double(img) / max(img(:));

% Choose the seed point as the center of the image
seed = [round(size(img, 1) / 2), round(size(img, 2) / 2)];

% Intensity difference threshold
threshold = 0.1;

% Calculate the binary mask based on the intensity difference threshold
mask = abs(normalized_img - normalized_img(seed(1), seed(2))) <= threshold;

% Fill holes in the binary mask
mask = imfill(mask, 'holes');

% Run region growing on the binary mask
segmented_rg = region_growing(mask, seed, threshold);


end



function segmented = region_growing(img, seed, threshold)
    [rows, cols] = size(img);
    segmented = zeros(rows, cols);
    visited = false(rows, cols);
    queue = [seed];
    segmented(seed(1), seed(2)) = img(seed(1), seed(2));
    visited(seed(1), seed(2)) = true;

    while ~isempty(queue)
        current = queue(1, :);
        queue(1, :) = [];
        
        neighbors = [
            current(1) - 1, current(2);
            current(1) + 1, current(2);
            current(1), current(2) - 1;
            current(1), current(2) + 1;
        ];
        
        for i = 1:size(neighbors, 1)
            r = neighbors(i, 1);
            c = neighbors(i, 2);
            
            if r >= 1 && r <= rows && c >= 1 && c <= cols && ~visited(r, c)
                visited(r, c) = true;
                intensity_diff = abs(double(img(r, c)) - double(img(current(1), current(2))));
                
                if intensity_diff <= threshold
                    segmented(r, c) = img(r, c);
                    queue = [queue; [r, c]];
                end
            end
        end
    end
end
