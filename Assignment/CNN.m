% Load the MRI data from the Brain.mat file
load('Brain.mat');
MRI_data = double(T1);

% Normalize the MRI data to have zero mean and unit variance
mean_val = mean(MRI_data(:));
std_val = std(MRI_data(:));
MRI_data = (MRI_data - mean_val) / std_val;

% Reshape the MRI data into a 2D matrix and shuffle the order
X = reshape(MRI_data, [], size(MRI_data,3));
idx = randperm(size(X,2));
X = X(:, idx);

% Create the corresponding label matrix for the slices
% Here, we assume that the labels are stored in a separate matrix called 'Labels'
labels = reshape(label, [], size(label,3));
Y = labels(:, idx);

% Split the data into training, validation, and testing sets
num_train = 0.7*size(X,2);
num_valid = 0.1*size(X,2);
num_test = size(X,2) - num_train - num_valid;
X_train = X(:, 1:num_train);
Y_train = Y(:, 1:num_train);
X_valid = X(:, num_train+1:num_train+num_valid);
Y_valid = Y(:, num_train+1:num_train+num_valid);
X_test = X(:, num_train+num_valid+1:end);
Y_test = Y(:, num_train+num_valid+1:end);

% Import the required functions and classes
import matlab.net.*
import matlab.net.http.*
import matlab.net.http.field.*

% Load and prepare the data, as previously described
% ... (previous code snippet)

% Define the CNN architecture
layers = [
    imageInputLayer([size(MRI_data,1) size(MRI_data,2) 1]) % assuming grayscale images

    convolution2dLayer(3, 32, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2, 'Stride', 2)
    
    convolution2dLayer(3, 64, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2, 'Stride', 2)
    
    convolution2dLayer(3, 128, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    
    transposedConv2dLayer(4, 64, 'Cropping', 'same', 'Stride', 2)
    batchNormalizationLayer
    reluLayer
    
    transposedConv2dLayer(4, 32, 'Cropping', 'same', 'Stride', 2)
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer(1, 5, 'Padding', 'same') % 5 output channels for the 5 distinct layers
    softmaxLayer
    pixelClassificationLayer
];

% Specify training options
options = trainingOptions('adam', ...
    'InitialLearnRate', 1e-3, ...
    'MaxEpochs', 100, ...
    'MiniBatchSize', 16, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {X_valid, Y_valid}, ...
    'ValidationFrequency', 30, ...
    'Verbose', false, ...
    'Plots', 'training-progress');

% Train the CNN
net = trainNetwork(categorical(T1), categorical(label), layers, options);

% Evaluate the trained CNN on the test set
Y_pred = classify(net, X_test);
accuracy = sum(Y_pred == Y_test) / numel(Y_test);
disp(['Test Accuracy: ' num2str(accuracy)]);

