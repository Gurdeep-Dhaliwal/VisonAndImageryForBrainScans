load ('Brain.mat'); % Load the MRI Data

num_clusters = 1;
accuracy_table = zeros(2, 10);
percentages = [0, 6, 14, 31, 50.5, 65];

for i = 1:10
    slice = T1(:,:,i);
    label_slice = label(:,:,i);

    % Preprocessing
    img = imadjust(slice, stretchlim(slice), []);

    % Thresholding
    threshold = graythresh(img);
    bw = imbinarize(img, threshold);

    % Morphological operations
    se = strel('disk', 5);
    bw = imopen(bw, se);
    bw = imfill(bw, 'holes');
    bw = bwareaopen(bw, 50);

    % Labeling
    cc = bwconncomp(bw);
    stats = regionprops(cc, 'Area');
    [~, idx] = sort([stats.Area], 'descend');
    bw = ismember(labelmatrix(cc), idx(1:num_clusters));

    % Accuracy calculation
    layers_acc = dice(double(bw), double(label_slice));
    accuracy_table(1, i) = mean(layers_acc);
end

disp(accuracy_table);