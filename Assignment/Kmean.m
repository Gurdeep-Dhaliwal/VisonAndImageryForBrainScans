function Kmean2D = Kmean(img, num_clusters)
    % Create matrices for k-means clustering and initialize all elements to zero
    layersK = zeros(size(img));
    layersKmean = layersK;
    
    % Set the number of clusters to the input value of num_clusters
    k = num_clusters;
    
    % Reshape input image into a column vector
    X = reshape(img, [], 1);
    
    % Identify background pixels (with intensity value of 0)
    bg_pixels = (X == 0);
    
    % Perform k-means clustering on the column vector X using k clusters
    [idx, C] = kmeans(X, k, 'Distance', 'sqEuclidean', 'Replicates', 5);
    
    % Sort cluster centers in ascending order
    [~, I] = sort(C);
    
    % Assign layer values to each pixel based on its cluster index
    for j = 1:k
        layersKmean(idx == I(j)) = j;
    end
    
    % Return k-means clustered image
    Kmean2D = layersKmean;
end
