load('Brain.mat');

function segmented_image = kmeans_segmentation(image, num_clusters)
    reshaped_image = double(reshape(image, [], 1));
    [cluster_idx, ~] = kmeans(reshaped_image, num_clusters, 'Distance', 'sqEuclidean', 'Replicates', 5);
    segmented_image = reshape(cluster_idx, size(image));
end


function dsc = dice_similarity_coefficient(A, B)
    intersection = sum(sum(A & B));
    union = sum(sum(A | B));
    dsc = 2 * intersection / (sum(sum(A)) + sum(sum(B)));
end


num_slices = 10;
num_clusters = 6; % Including air as a label
segmented_images = cell(1, num_slices);

for i = 1:num_slices
    segmented_images{i} = kmeans_segmentation(T1(:,:,i), num_clusters);
end

total_dsc = 0;
num_labels = 6;

for i = 1:num_slices
    for label = 0:(num_labels-1)
        A = segmented_images{i} == label;
        B = PreSegmented(:,:,i) == label;
        dsc = dice_similarity_coefficient(A, B);
        total_dsc = total_dsc + dsc;
    end
end

average_dsc = total_dsc / (num_slices * num_labels);
fprintf('Average DSC for K-means segmentation: %.4f\n', average_dsc);







%Thresholding
    min_intensity = min(img(:));
    max_intensity = max(img(:));
    thresholds = min_intensity + (percentages / 100) * (max_intensity - min_intensity);
    layers = zeros(size(img));
    for j = 1:5
        layers(img > thresholds(j) & img <= thresholds(j+1)) = j;
    end
    %Thresholding accuracy
    layers_acc = dice(layers, double(label_slice));
    accuracy_table(2, i) = mean(layers_acc);





    % Preprocessing for the watershed method
    img3 = slice;
    img3 = imgaussfilt(img3, 2);
    img3 = imsharpen(img3, 'Amount', 0.5, 'Radius', 1.5);

    % Watershed segmentation
    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(img3), hy, 'replicate');
    Ix = imfilter(double(img3), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);

    
    % Marker computation
    se = strel('disk', 20);
    Io = imopen(img3, se);
    Ie = imerode(img3, se);
    Iobr = imreconstruct(Ie, img3);
    Ioc = imclose(Io, se);
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    fgm = imregionalmax(Iobrcbr);
    I2 = imcomplement(img3);
    fgm4 = bwareaopen(fgm, 20);
    I3 = imimposemin(I2, fgm4 | fgm);
    L = watershed(I3);

    imshow(L,[])
    






















    load ('Brain.mat'); % Load the MRI Data

%T1 weighted MRI data at 10 consecutive slices
%{
for i = 1:10                                     %Looping through every slice
    figure(1);
    subplot(2,5,i);                             %Creating a subplot to display all slices
    slice = T1(:,:,i);                          %Extracting the slices
    imagesc(slice);                           %Display the slices
    colormap(gray);                         %Using the grayscale colormap
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;                           %Setting the axis to have equal ratios and removing all axis lables
    title(sprintf("Slice %d", i));       %Title for the subplot
end
sgtitle("All Slices");                      %Title for the figure to display brain slices

%Pre segmented data
for i = 1:10                                     %Looping through every pre-segemented slice
    figure(2);
    subplot(2,5,i);                             %Creating a subplot to display all pre-segemented slices
    slice = label(:,:,i);                       %Extracting the pre-segemented slices
    imagesc(slice);                           %Display the pre-segemented slices
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;                           %Setting the axis to have equal ratios and removing all axis lables
    title(sprintf("Slice %d", i));       %Title for the subplot
end
sgtitle("Pre-segemented Slices");  %Title for the figure to display pre-segemented slices
%}



num_clusters = 5;
accuracy_table = zeros(4, 10);
percentages = [0, 7, 19, 36, 58.4, 74.1];


% Background removal threshold
bg_threshold = 31000;


for i = 1:10
    slice = T1(:,:,i);
    label_slice = label(:,:,i);

    slice(slice < bg_threshold) = 0;

    %figure(2);
    %g = histogram(slice);
  
    %Pre-Processing
    %img = imsharpen(slice);
    %img = medfilt2(img);

    img = medfilt2(slice);
    img = imgaussfilt(img, 2);
    img = imsharpen(img, 'Amount', 0.7, 'Radius', 2);



    %Canny
    Canny = edge(img, 'canny');
    %Canny accuracy
    canny_acc = dice(Canny, logical(label_slice));
   
    %Manual Thresholding
    min_intensity = min(img(:));
    max_intensity = max(img(:));
    thresholds = min_intensity + (percentages / 100) * (max_intensity - min_intensity);
    layers3 = zeros(size(img));
    for j = 1:5
        layers3(img > thresholds(j) & img <= thresholds(j+1)) = j;
    end
    %Thresholding accuracy
    layers_acc2 = dice(layers3, double(label_slice));
    accuracy_table(1, i) = mean(layers_acc2);


    %Automatic Thresholding using Otsu
    thresholds = multithresh(slice, 5);                                                        
    thresholds = [0 , thresholds];
    layers = zeros(size(slice));
    for j = 1:5
        layers(slice > thresholds(j) & slice <= thresholds(j+1)) = j;
    end

    %Otsu Thresholding accuracy
    layers_acc = dice(layers, double(label_slice));
    accuracy_table(2, i) = mean(layers_acc);
    

    %Automatic Thresholding
    T = multithresh(img, 5);
    T = [0, T];
    layers2= zeros(size(img));
    for k = 1:5
        layers2(img > T(k) & img <= T(k+1)) = k;
    end

    layers2_acc = dice(layers2, double(label_slice));
    
   

    % K-means clustering
    layersK = zeros(size(img));
    layersKmean = layersK;
    k = num_clusters;
    X = reshape(img, [], 1);
    bg_pixels = (X == 0);
    [idx, C] = kmeans(X, k);
    [~, I] = sort(C);
    for j = 1:k
        layersKmean(idx == I(j)) = j;
    end

    layersKmean(bg_pixels) = 0;

    layersKmean_acc = dice(layersKmean, double(label_slice));
    accuracy_table(3,i) = mean(layersKmean_acc);


    % Watershed Algorithm
    img_d = double(img);
    img_grad = imgradient(img_d);
    img_grad_min = imimposemin(img_grad, Canny);
    layers_watershed = watershed(img_grad_min);
    layers_watershed(layers_watershed == 1) = 0;

    imshow(layers_watershed,[])

end

    
    %Displaying Original slices next to segemented slices
    
    %{
    %Displaying original slices
    subplot(4,10,(i-1)*4+1);                                  
    imagesc(slice);                                 
    axis image off;                                 
    title(sprintf('Slice %d - Original', i)); 

    subplot(4,10,(i-1)*4+2);                                  
    imagesc(Canny);                                 
    axis image off;
    title(sprintf('Slice %d - Canny', i))

    subplot(4,10,(i-1)*4+3);                                  
    imagesc(layers);                                 
    axis image off;
    title(sprintf('Slice %d - Threshold', i))

    subplot(4,10,(i-1)*4+4);                                  
    imagesc(segmented_image);                                 
    axis image off;
    title(sprintf('Slice %d - Kmean', i))
    %}


disp(accuracy_table);
writematrix(accuracy_table, 'accuracy_table.csv');

%{
%Displaying Original slices next to segemented slices

    %Displaying original slices
    subplot(4,1,1);                                  
    imagesc(slice);        
    axis image off;                                 
    title(sprintf('Slice %d - Original', i)); 

    subplot(4,1,2);                                  
    imagesc(Canny); 
    axis image off;
    title(sprintf('Slice %d - Canny', i))

    subplot(4,1,3);                                  
    imagesc(layers); 
    axis image off;
    title(sprintf('Slice %d - Threshold', i))

    subplot(4,1,4);                                  
    imagesc(layersKmean);       
    axis image off;
    title(sprintf('Slice %d - Kmean', i))
%}    



















load ('Brain.mat'); % Load the MRI Data

%T1 weighted MRI data at 10 consecutive slices
%{
for i = 1:10                                     %Looping through every slice
    figure(1);
    subplot(2,5,i);                             %Creating a subplot to display all slices
    slice = T1(:,:,i);                          %Extracting the slices
    imagesc(slice);                           %Display the slices
    colormap(gray);                         %Using the grayscale colormap
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;                           %Setting the axis to have equal ratios and removing all axis lables
    title(sprintf("Slice %d", i));       %Title for the subplot
end
sgtitle("All Slices");                      %Title for the figure to display brain slices

%Pre segmented data
for i = 1:10                                     %Looping through every pre-segemented slice
    figure(2);
    subplot(2,5,i);                             %Creating a subplot to display all pre-segemented slices
    slice = label(:,:,i);                       %Extracting the pre-segemented slices
    imagesc(slice);                           %Display the pre-segemented slices
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;                           %Setting the axis to have equal ratios and removing all axis lables
    title(sprintf("Slice %d", i));       %Title for the subplot
end
sgtitle("Pre-segemented Slices");  %Title for the figure to display pre-segemented slices
%}



num_clusters = 5;
accuracy_table = zeros(4, 10);
percentages = [0, 3, 18, 34.7, 54.3, 69.5];

% Background removal threshold
bg_threshold = 31000;


for i = 1:10
    slice = T1(:,:,i);
    label_slice = label(:,:,i);

    slice(slice < bg_threshold) = 0;

    %figure(2);
    %g = histogram(slice);
  
    %Pre-Processing
    img = imsharpen(slice);
    img = medfilt2(img);

    img2 = slice;
    img2 = imgaussfilt(img2, 2);
    img2 = imsharpen(img2, 'Amount', 0.5, 'Radius', 1.5);

    %Canny
    Canny = edge(img, 'canny');
    %Canny accuracy
    canny_acc = dice(Canny, logical(label_slice));
   
    %Manual Thresholding
    min_intensity = min(img(:));
    max_intensity = max(img(:));
    thresholds = min_intensity + (percentages / 100) * (max_intensity - min_intensity);
    layers3 = zeros(size(img));
    for j = 1:5
        layers3(img > thresholds(j) & img <= thresholds(j+1)) = j;
    end
    %Thresholding accuracy
    layers_acc2 = dice(layers3, double(label_slice));
    accuracy_table(1, i) = mean(layers_acc2);


    %Automatic Thresholding using Otsu
    thresholds = multithresh(slice, 5);                                                        
    thresholds = [0 , thresholds];
    layers = zeros(size(slice));
    for j = 1:5
        layers(slice > thresholds(j) & slice <= thresholds(j+1)) = j;
    end

    %Thresholding accuracy
    layers_acc = dice(layers, double(label_slice));
    accuracy_table(2, i) = mean(layers_acc);
    

    %Automatic Thresholding
    T = multithresh(img2, 5);
    T = [0, T];
    layers2= zeros(size(img2));
    for k = 1:5
        layers2(img2 > T(k) & img2 <= T(k+1)) = k;
    end

    layers2_acc = dice(layers2, double(label_slice));
    
   

    % K-means clustering
    layersK = layers;
    layersKmean = layersK;
    k = num_clusters;
    X = reshape(slice, [], 1);
    bg_pixels = (X == 0);
    [idx, C] = kmeans(X, k);
    [~, I] = sort(C);
    for j = 1:k
        layersKmean(idx == I(j)) = j;
    end

    layersKmean(bg_pixels) = 0;

    layersKmean_acc = dice(layersKmean, double(label_slice));
    accuracy_table(3,i) = mean(layersKmean_acc);



    % Preprocessing for the watershed method
    img3 = slice;
    img3 = imgaussfilt(img3, 2);
    img3 = imsharpen(img3, 'Amount', 0.5, 'Radius', 1.5);

    % Watershed segmentation
    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(img3), hy, 'replicate');
    Ix = imfilter(double(img3), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);

    
    % Marker computation
    se = strel('disk', 20);
    Io = imopen(img3, se);
    Ie = imerode(img3, se);
    Iobr = imreconstruct(Ie, img3);
    Ioc = imclose(Io, se);
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    fgm = imregionalmax(Iobrcbr);
    I2 = imcomplement(img3);
    fgm4 = bwareaopen(fgm, 20);
    I3 = imimposemin(I2, fgm4 | fgm);
    L = watershed(I3);

    imshow(I3,[])
    


    % Watershed segmentation accuracy
    watershed_acc = dice(watershed_labels, double(label_slice));
    accuracy_table(4, i) = mean(watershed_acc);


end

