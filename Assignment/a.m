%Displaying original slices
    %subplot(10,3,(i-1)*3+1);                                  % Create a subplot for the original slice
    %imagesc(slice);                                 % Display the slice
    %axis image off;                                  % Set the axes to have equal ratios
    %title(sprintf('Slice %d - Original', i)); % Add title to the subplot


    %Create subplot for thresholded slices
    %subplot(10,3,(i-1)*3+2);                                   
    %imagesc(layers);                       
    %colormap(gray); 
    %axis image off;                                    
    %title(sprintf('Slice %d - layers', i)); 

    
    %subplot(10,3,(i-1)*3+3);                                   
    %imagesc(T1_Canny);                       
    %colormap(gray); 
    %axis image off;                                    
    %title(sprintf('Slice %d - img_clean', i)); 






      %
    gmag = imgradient(img);
    gmag = imsharpen(gmag);
    gmag = medfilt2(gmag);


    WSmin_intensity = min(gmag(:));
    WSmax_intensity = max(gmag(:));
    WSthresholds = WSmin_intensity + (WSpercentages / 100) * (WSmax_intensity - WSmin_intensity);
    WSlayers = zeros(size(gmag));
    for k = 1:5
        WSlayers(gmag > WSthresholds(k) & gmag <= WSthresholds(k+1)) = k;
    end

    img_clean = imfill(WSlayers, 'holes');

    clean_acc = dice(img_clean, double(label_slice));
    accuracy_table(2,i) = mean(clean_acc);
