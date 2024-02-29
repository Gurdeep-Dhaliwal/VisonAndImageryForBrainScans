load ('Brain.mat'); % Load the MRI Data

%% Displaying T1- Data to be segmented and  labels- Pre-segmented Data
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



%% Setting values
num_clusters = 5;
accuracy_table2D = zeros(5, 10);
accuracy_table3D = zeros(5,10);
percentages = [0, 7, 17.6, 37, 57.4, 74.3];
percentages3D = [0, 5, 16, 32.4, 52.6, 67];

% Background removal threshold
bg_threshold = 31000;

for i = 1:10                                     %Looping through the 10 slices
%%    2D segmenation
    slice = T1(:,:,i);                          %Extrating the slices
    label_slice = label(:,:,i);             %Extracting pre-segmented slices

    slice(slice < bg_threshold) = 0;  %Setting background threshold
    bg_pixels = (slice == 0);             

    %figure(2);
    %g = histogram(slice);
  
    %Pre-Processing
    img = medfilt2(slice,[6,5]);         %To remove noise and irregularites in the image
    img = imgaussfilt(img, 0.1);        %Gaussian fitler to smooth image
    img = imsharpen(img, 'Amount', 0.01, 'Radius', 0.9);  %Sharpen image
    img = double(img);

%%    2D Canny
    Canny = edge(img, 'canny');
    %Canny accuracy
    canny_acc = dice(Canny, logical(label_slice));

 %%  Watershed segmentation
    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(img), hy, 'replicate');
    Ix = imfilter(double(img), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);

    % Marker computation
    se = strel('disk', 20);
    Io = imopen(img, se);
    Ie = imerode(img, se);
    Iobr = imreconstruct(Ie, img);
    Ioc = imclose(Io, se);
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    fgm = imregionalmax(Iobrcbr);
    I2 = imcomplement(img);
    fgm4 = bwareaopen(fgm, 20);
    I3 = imimposemin(I2, fgm4 | fgm);
    L = watershed(I3);
  
%%  2D  Manual Thresholding
    [MT, thresholds] = ManualThresholding(img, percentages);
    Mlayers = MT;

    Mlayers(bg_pixels) =  0;

    %Thresholding accuracy
    layers_acc2 = dice(Mlayers, double(label_slice));
    accuracy_table2D(1, i) = mean(layers_acc2(2:end));

%%    2D Automatic Thresholding using Otsu
    Alayers = AutomaticThresholding(img);
    Alayers(bg_pixels) = 0;

    %Otsu Thresholding accuracy
    layers_acc = dice(Alayers, double(label_slice));
    accuracy_table2D(2, i) = mean(layers_acc(2:end));
   
%%    2D K-means clustering
    layersKmean = Kmean(img, num_clusters);
    layersKmean(bg_pixels) = 0;

    layersKmean_acc = dice(layersKmean, double(label_slice));
    accuracy_table2D(3,i) = mean(layersKmean_acc(2:end));

%% 2D Region Growing
    % Using thresholding as an initial segmentation
    init_seg = zeros(size(img));
    for j = 1:5
        init_seg(slice > thresholds(j) & slice <= thresholds(j+1)) = j;
    end

    seed = [round(size(img, 1) / 2), round(size(img, 2) / 2)];
    segmented_rg = RegionGrowing(slice, init_seg, seed, 0.17);
    segmented_rg(bg_pixels) = 0;

    % Region Growing accuracy
    segmented_rg_acc = dice(segmented_rg, double(label_slice));
    accuracy_table2D(4, i) = mean(segmented_rg_acc(2:end));

    % Use K-means segmentation as an initial segmentation for region growing
    init_seg_kmeans = layersKmean;

    seed = [round(size(img, 1) / 2), round(size(img, 2) / 2)];
    segmented_rg_kmeans = RegionGrowing(slice, init_seg_kmeans, seed, 0.15);
    segmented_rg_kmeans(bg_pixels) = 0;

    % Combined K-means and Region Growing accuracy
    segmented_rg_kmeans_acc = dice(segmented_rg_kmeans, double(label_slice));
                             accuracy_table2D(5, i) = mean(segmented_rg_kmeans_acc(2:end));


%%   3D segmentation

    %Pre-Processing
    T1(T1<bg_threshold) = 0;
    bg_pixels3D = (T1 == 0);
    img3D = medfilt3(T1,[5,3,1]);
    img3D = double(imgaussfilt3(img3D, 0.01));
    img3DD = double(T1);
    
    

%%    3D Manual Thresholding
    %Thresholding accuracy
    [MT, thresholds] = ManualThresholding(T1, percentages3D);
    M3Dlayers = MT;
    M3Dthresholds = thresholds;
    M3Dlayers(bg_pixels3D) =  0;

    %Thresholding accuracy
    M3Dlayers_acc = dice(M3Dlayers, double(label));
    accuracy_table3D(1, i) = mean(M3Dlayers_acc(2:end));

%%    3D Automatic Thresholding using Otsu
    A3Dlayers = AutomaticThresholding(T1);

    A3Dlayers_acc = dice(A3Dlayers, double(label));
    accuracy_table3D(2,i) = mean(A3Dlayers_acc(2:end));



%% 3D K-means clustering
    layersK3Dmean = Kmean3D(img3D, num_clusters);

    layersK3Dmean(bg_pixels3D) = 0;

    layersK3Dmean_acc = dice(layersK3Dmean, double(label));
    accuracy_table3D(3,i) = mean(layersK3Dmean_acc(2:end));
  

%% 3D K-means clustering followed by 3D Region Growing
    %Using Thresholding
    init_seg_threshold_3D = M3Dlayers;

    seed_3d = [size(img3D, 1) / 2, size(img3D, 2) / 2, size(img3D, 3) / 2];
    segmented_rg_threshold_3D = region_growing_3d(T1, init_seg_threshold_3D, seed_3d, 1);

    segmented_rg_threshold_3D_acc = dice(segmented_rg_threshold_3D, double(label));
    accuracy_table3D(4, i) = mean(segmented_rg_threshold_3D_acc(2:end));

    % Use K-means segmentation as an initial segmentation for region growing
    init_seg_kmeans_3d = layersK3Dmean;

    segmented_rg_kmeans_3d = region_growing_3d(T1, init_seg_kmeans_3d, seed_3d, 0.15);

    % Combined K-means and Region Growing accuracy
    segmented_rg_kmeans_3d_acc = dice(segmented_rg_kmeans_3d, double(label));
    accuracy_table3D(5, i) = mean(segmented_rg_kmeans_3d_acc(2:end));

    
%%    Displaying Original slices next to segemented slices(2D)
    
    %{
    %Displaying original slices
    subplot(6,10,(i-1)*6+1);                                  
    imagesc(slice);                                 
    axis image off;                                 
    title(sprintf('Slice %d - Original', i)); 

    subplot(6,10,(i-1)*6+2);                                  
    imagesc(Canny);                                 
    axis image off;
    title(sprintf('Slice %d - Canny', i))

    subplot(6,10,(i-1)*6+3);                                  
    imagesc(Mlayers);                                 
    axis image off;
    title(sprintf('Slice %d - Manual Threshold', i))

    subplot(6,10,(i-1)*6+4);                                  
    imagesc(Alayers);                                 
    axis image off;
    title(sprintf('Slice %d - Automatic Threshold', i))

    subplot(6,10,(i-1)*6+5);                                  
    imagesc(layersKmean);                                 
    axis image off;
    title(sprintf('Slice %d - Kmean', i))

    subplot(6,10,(i-1)*6+6);                                  
    imagesc(segmented_rg);                                 
    axis image off;
    title(sprintf('Slice %d - Reigon Growing', i))
    %}

    %% Displaying Original slice next to segemented slice(2D, 1st slice)
    if i == 1
        subplot(6,10,1);                                  
    imagesc(slice);
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;                                 
    title(sprintf('Slice %d - Original', 1)); 

    subplot(6,10,2);                                  
    imagesc(Mlayers); 
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Manual Threshold', 1))

    subplot(6,10,3);                                  
    imagesc(Alayers);
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Automatic Threshold', 1))

    subplot(6,10,4);                                  
    imagesc(layersKmean);
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Kmean', 1))

    subplot(6,10,5);                                  
    imagesc(segmented_rg); 
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Reigon Growing T', 1))
   
    subplot(6,10,6);                                  
    imagesc(segmented_rg_kmeans); 
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Reigon Growing KM', 1))
    end
 end

%% Displaying Accuracy Tables
disp("Title : Accuracy Table 2D")
disp(accuracy_table2D);
writematrix(accuracy_table2D, 'accuracy_table2D.csv');

disp("Title: Accuracy Table 3D")
disp(accuracy_table3D);
writematrix(accuracy_table3D, 'accuracy_table3D.csv');

%% Displaying Original slice next to segemented slice(3D, 1st slice)
%{
    %Displaying original slices
    subplot(6,10,1);                                  
    imagesc(T1(:,:,1));
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;                                 
    title(sprintf('Slice %d - Original', 1)); 

    subplot(6,10,2);                                  
    imagesc(M3Dlayers(:,:,1)); 
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Manual Threshold', 1))

    subplot(6,10,3);                                  
    imagesc(A3Dlayers(:,:,1));
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Automatic Threshold', 1))

    subplot(6,10,4);                                  
    imagesc(layersK3Dmean(:,:,1));
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Kmean', 1))

    subplot(6,10,5);                                  
    imagesc(segmented_rg_threshold_3D(:,:,1)); 
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Reigon Growing T', 1))
   
    subplot(6,10,6);                                  
    imagesc(segmented_rg_kmeans_3d(:,:,1)); 
    c = colorbar('southoutside');      % Add a colorbar to show intensity
    c.Label.String = 'Tissue Layers';% Labelling the colorbar
    axis image off;
    title(sprintf('Slice %d - Reigon Growing KM', 1))
%}