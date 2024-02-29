function segmented = RegionGrowing(img, init_seg, seed, threshold)
    % Get the number of rows and columns in the input image
    [rows, cols] = size(img);
    
    % Initialize the output segmented mask to be the same as the input init_seg
    segmented = init_seg;
    
    % Create a boolean matrix to keep track of which pixels have been visited
    visited = false(rows, cols);
    
    % Create a queue data structure containing the seed point and its initial segmentation value
    queue = [seed, init_seg(seed(1), seed(2))];
    
    % Mark the seed point as visited
    visited(seed(1), seed(2)) = true;

    % Loop until the queue is empty
    while ~isempty(queue)
        % Remove the first element from the queue and assign it to the current variable
        current = queue(1, :);
        queue(1, :) = [];
        
        % Compute the four neighboring pixels of the current pixel
        neighbors = [            
            current(1) - 1, current(2);            
            current(1) + 1, current(2);            
            current(1), current(2) - 1;          
            current(1), current(2) + 1;       
            ];
        
        % Loop through each neighboring pixel
        for i = 1:size(neighbors, 1)
            r = neighbors(i, 1);
            c = neighbors(i, 2);
            
            % Check whether the neighboring pixel is within the bounds of the input image and whether it has already been visited
            if r >= 1 && r <= rows && c >= 1 && c <= cols && ~visited(r, c)
                % Mark the neighboring pixel as visited
                visited(r, c) = true;
                
                % Compute the intensity difference between the neighboring pixel and the current pixel
                intensity_diff = abs(double(img(r, c)) - double(img(current(1), current(2))));
                
                % Check whether the intensity difference is less than or equal to the input threshold
                if intensity_diff <= threshold
                    % Assign the neighboring pixel the same segmentation value as the current pixel and add it to the queue
                    segmented(r, c) = current(3);
                    queue = [queue; [r, c, current(3)]];
                end
            end
        end
    end
    
    % Perform median filtering on the segmented image to remove any noise
    segmented = medfilt2(segmented, [5, 5]);
end

