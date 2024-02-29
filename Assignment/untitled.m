load ('Brain.mat'); % Load the MRI Data

%Displaying the data

%T1 weighted MRI data at 10 consecutive slices
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

%Applying Segmentation algorithms
num_clusters = 5;                           %Defining the number of k-means clusters
figure;                                             % Create a new figure window
accuracy_table = zeros(5, 10);       % Initialize accuracy table

for i = 1:10
    slice = T1(:,:,i);
    label_slice = label(:,:,i);

    %Applying k-mean clustering algorithm
    [xlabel, ylabel] = meshgrid(1:size(slice,2), 1:size(slice,1));
    data = [xlabel(:) ylabel(:) double(slice(:))];
    [~, centroids] = kmeans(data, num_clusters);
    [~, cluster_labels] = pdist2(centroids, data , 'euclidean', 'Smallest',1);
    T1_Kmeans = reshape(cluster_labels, size(slice));
    layer_labels = zeros(size(T1_Kmeans));
    for j = 1:num_clusters
        layer_labels(T1_Kmeans == j) = j;
    end

    %Appyling Canny edge detection 
    T1_Canny = edge(slice, 'canny');


    %Appyling Otsu's thresholding algorithm
    threshold = graythresh(slice);
    T1_Otsu = imbinarize(slice, threshold);

    % Applying thresholding
    if slice>0

    end

    %k-means accuracy
    kmeans_acc = dice(T1_Kmeans, double(label_slice));
    accuracy_table(2, i) = mean(kmeans_acc);

    %Canny accuracy
    canny_acc = dice(T1_Canny, logical(label_slice));
    accuracy_table(3, i) = mean(canny_acc);

    %Otsu accuracy
    level = graythresh(slice);
    T1_otsu = imbinarize(slice, level);

    %Thresholding accuracy
    thresholding_acc = dice(T1_Thresholding, logical(label_slice));
    accuracy_table(5,i) = mean(thresholding_acc);


    histSlice = histeq(slice);
    threshold1 = 0.02; % Threshold value for air and skin/scalp
    threshold2 = 0.04; % Threshold value for skull
    air_skin = histSlice < threshold1; % Air and skin/scalp are below threshold1
    skull = (histSlice >= threshold1) & (histSlice < threshold2); % Skull is between threshold1 and threshold2
    histSlice(histSlice >= threshold2) = threshold2; % Set intensities above threshold2 to threshold2
    % Note: We set the intensities above threshold2 to threshold2 to avoid including the brain tissue in the skull region
    
    % Region-growing: Segment out the CSF, gray matter, and white matter regions using region-growing
    csf = regiongrowing(histSlice, 200, 300, 3); % Seed point at (200, 300), threshold set to 3
    gray = regiongrowing(histSlice, 150, 250, 4); % Seed point at (150, 250), threshold set to 4
    white = regiongrowing(histSlice, 250, 200, 4); % Seed point at (250, 200), threshold set to 4

    figure;
    subplot(2, 3, 1); imshow(air_skin); title('Air and Skin/Scalp');
    subplot(2, 3, 2); imshow(skull); title('Skull');
    subplot(2, 3, 3); imshow(csf); title('CSF');
    subplot(2, 3, 4); imshow(gray); title('Gray Matter');
    subplot(2, 3, 5); imshow(white); title('White Matter');
    

    %Displaying Original slices next to segemented slices

    %Displaying original slices
    subplot(10,5,(i-1)*5+1);                                  % Create a subplot for the original slice
    imagesc(slice);                                 % Display the slice
    axis image off;                                  % Set the axes to have equal ratios
    title(sprintf('Slice %d - Original', i)); % Add title to the subplot

    %Create subplot for k-means segmentation slices
    subplot(10,5,(i-1)*5+2);                                    % Create a subplot for the k-means segmentation
    imagesc(T1_Kmeans);                       % Display the segmented image
    colormap(lines(num_clusters)); 
    axis image off;                                    % Set the axes to have equal aspect ratio
    title(sprintf('Slice %d - k-means', i)); % Add a title to the subplot

    %Create subplot for Canny segmentation slices
    subplot(10,5,(i-1)*5+3);                                    % Create a subplot for the Canny segmentation
    imagesc(T1_Canny);                       % Display the segmented image
    axis image off;                                    % Set the axes to have equal aspect ratio
    title(sprintf('Slice %d - Canny', i)); % Add a title to the subplot

    %Create subplot for Globabl thresholding segmentation slices
    subplot(10,5,(i-1)*5+4);                                    % Create a subplot for the globabl thresholding segmentation
    imagesc(T1_Thresholding);                       % Display the segmented image
    axis image off;                                    % Set the axes to have equal aspect ratio
    title(sprintf('Slice %d - GT', i)); % Add a title to the subplot
end

%Analyze the data
 disp(accuracy_table);

