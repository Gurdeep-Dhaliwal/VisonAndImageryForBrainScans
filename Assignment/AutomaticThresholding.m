function AT = AutomaticThresholding(img)
    % Calculate 5 threshold values for the input image using multithresh
    ATthresholds = multithresh(img, 5);
    
    % Add 0 to the beginning of the threshold values to ensure that all pixels with intensities below the first threshold are included in the lowest threshold layer
    ATthresholds = [0, ATthresholds];
    
    % Create matrix for thresholded image and initialize all elements to zero
    ATlayers = zeros(size(img));
    
    % Loop through each threshold value
    for j = 1:5
        % Set elements in ATlayers that fall within the range (ATthresholds(j), ATthresholds(j+1)] to the value j
        ATlayers(img > ATthresholds(j) & img <= ATthresholds(j+1)) = j;
    end
    
    % Return thresholded image
    AT = ATlayers;
end
