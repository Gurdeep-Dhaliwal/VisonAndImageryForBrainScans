%Read Image
Fish_Vis = imread('fish-vis.tif');
Fish_CFP = imread('fish-cfp-1.tif');

%Covert to GreyScale
Fish_CFG = im2gray(Fish_CFP);

%Use cpselect() and save control points
fixedPoints = [250 200; 250 150; 250 100; 150 200; 150 150; 150 100];
movingPoints = [250 200; 250 150; 250 100; 150 200; 150 150; 150 100];
h = cpselect(Fish_CFP,Fish_Vis, movingPoints,fixedPoints,'Wait',false);
%[movingPoints,fixedPoints] = cpselect(Fish_CFG,Fish_Vis ,"Wait",true);

%Determine parameteres of transformation 
%tform = fitgeotform2d(movingPoints,fixedPoints,"similarity");
tform = fitgeotrans(movingPoints,fixedPoints,"projective");

%Transform the Image
registered_image = imwarp(Fish_CFP,tform,'FillValues',0,'OutputView', imref2d(size(Fish_Vis))); 

%Display both registered and base Image 
imshowpair(Fish_Vis,registered_image)

alpha=0.6; 
hold on 
h = imshow(Fish_Vis, gray(256)); 
set(h, 'AlphaData', alpha); 
