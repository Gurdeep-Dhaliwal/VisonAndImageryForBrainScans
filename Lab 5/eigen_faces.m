% This script will help you understand how face recognistion works by the
% use of linear dimensionality reduction through the use of PCA

%% First, let us read all the images and plot a few of them
% Note that raw_images is a cell variable with 86 images (2-4 images of
% multiple people
load raw_images.mat
figure('units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
hold on
k = 1;
for i = 1 : 7 : 70
    subplot(2,5,k)
    imagesc(raw_images{i});
    colormap gray; axis equal; axis off
    title(['Image Number = ' num2str(i)]);
    k = k+1;
end

%% Convert data to 2-D matrix, one row per image
% This means we now have a 2D array, each row an image (an observation) and
% each column being a data point (pixel value)
[x, y] = size(raw_images{1});
images = zeros(length(raw_images), x * y);
for i = 1:length(raw_images)
    images(i,:) = reshape(double(raw_images{i}), 1, x * y);
end


% save image resolution for later!
im_res = [x y];
clear x y raw_images i

%% Eigen Vectors (also known as eigenfaces
% next we will go through each step of calculating the PCA
% note the variable images, has 86 observations, with 3000 points each and
% to view, for example observation 11, we can say
figure('units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
hold on
imagesc(reshape(images(11,:),im_res));
set(gca,'Ydir','reverse');
axis equal; axis tight; colormap gray

% Typically you would subtract the mean of each observation, then form the
% covariance matrix and then caluclate the Principal components, but our
% data space is large 3000 points so a big computational problem.
% We would normally do the following, but for higher resolution we need to
% think differently! Can you think how?
% Step 1: Calculate co-variance
x = images;
mean_x = mean(x,2);
r = x-repmat(mean_x,1,size(x,2));
conv_x = r'*r;

% GEt Eigenvalues and Eigenvectors 
[V,D] = eig(conv_x);  
Eig = diag(D);
% and sort them
[~,idx] = sort(Eig,'descend') ;
values = Eig(idx);
vectors = V(:,idx);
clear V D idx

% plot the first 6 PCs
figure('units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
hold on
k = 1;
for i  = 1 : 6
    subplot(2,3,k)
    imagesc(reshape(vectors(:,i),im_res));axis equal;
    axis tight; colormap gray; axis off
    title(['PC = ' num2str(i)]);
    k = k + 1;
end

%%
% each row is the weight (reduced dimension) of each original image
% so for each face, we can calculate the weight w.
w = (images)*vectors;

figure('units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
hold on
imagesc(reshape(vectors*w(61,:)',im_res));
axis equal; axis tight; colormap gray
title('Face 61, Using ALL of the eigenVectors');

% Now recover original data but only use for example 50 Dimensions, that is
% the first 50 EigenVectors
w1 = w; 
w1(:,30:end) = 0;
figure('units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
hold on
imagesc(reshape(vectors*w1(61,:)',im_res));axis equal; axis tight; colormap gray
title('Face 61, Using only the first 30 eigenVectors');

% and the original
figure('units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
hold on
imagesc(reshape(images(61,:),im_res));axis equal; axis tight; colormap gray
title('Face 61, original');

