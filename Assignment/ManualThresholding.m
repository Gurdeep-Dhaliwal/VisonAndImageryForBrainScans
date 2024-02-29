function [MT, thresholds] = ManualThresholding(img, percentages)
    % Determine minimum and maximum intensity values in input image
    min_intensity = min(img(:));
    max_intensity = max(img(:));
    
    % Calculate threshold values using input percentages
    thresholds = min_intensity + (percentages / 100) * (max_intensity - min_intensity);
    
    % Create matrix for thresholded image and initialize all elements to zero
    Mlayers = zeros(size(img));
    
    % Loop through each threshold value
    for j = 1:5
        % Set elements in Mlayers that fall within the range (threshold(j), threshold(j+1)] to the value j
        Mlayers(img > thresholds(j) & img <= thresholds(j+1)) = j;
    end
    
    % Return thresholded image and threshold values
    MT = Mlayers;
end
