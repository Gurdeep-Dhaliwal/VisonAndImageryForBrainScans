function Kmean3D = Kmean3D(img, num_clusters)
    % Create matrices for k-means clustering and initialize all elements to zero
    layersK3D = zeros(size(img));
    layersK3Dmean = layersK3D;
    
    % Set the number of clusters to the input value of num_clusters
    k3D = num_clusters;
    
    % Reshape input 3D image into a column vector
    X3D = reshape(img, [], 1);
    
    % Identify background pixels (with intensity value of 0)
    bg_pixels3D = (X3D == 0);
    
    % Perform k-means clustering on the column vector X using k clusters
    [idx, C] = kmeans(X3D, k3D,'Distance', 'sqEuclidean', 'Replicates', 5);
    
    % Sort cluster centers in ascending order
    [~, I] = sort(C);
    
    % Assign layer values to each pixel based on its cluster index
    for j = 1:k3D
        layersK3Dmean(idx == I(j)) = j;
    end
    
    % Reshape clustered image into its original 3D shape
    Kmean3D = reshape(layersK3Dmean, size(img));
end
