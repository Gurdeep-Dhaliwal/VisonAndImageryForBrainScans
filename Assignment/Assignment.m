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


figure;
num_clusters = 5;
accuracy_table = zeros(2, 10);

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

    %k-means accuracy
    kmeans_acc = dice(T1_Kmeans, double(label_slice));
    accuracy_table(2, i) = mean(kmeans_acc);



 %Displaying Original slices next to segemented slices

    %Displaying original slices
    subplot(2,10,(i-1)*2+1);                                  % Create a subplot for the original slice
    imagesc(slice);                                 % Display the slice
    axis image off;                                  % Set the axes to have equal ratios
    title(sprintf('Slice %d - Original', i)); % Add title to the subplot

    %Create subplot for k-means segmentation slices
    subplot(2,10,(i-1)*2+2);                                    % Create a subplot for the k-means segmentation
    imagesc(T1_Kmeans);                       % Display the segmented image
    colormap(lines(num_clusters)); 
    axis image off;                                    % Set the axes to have equal aspect ratio
    title(sprintf('Slice %d - k-means', i)); % Add a title to the subplot
end

disp(accuracy_table);
