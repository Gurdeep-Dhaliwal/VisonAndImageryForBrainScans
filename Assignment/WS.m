load("Brain.mat");

I = T1(:,:,1);
percentages = [0, 6, 15.5, 31, 53, 68];

gmag = imgradient(I);
imshow(gmag,[])
title("gmag")

L = watershed(gmag);
Lrgb = label2rgb(L);


se = strel("disk", 20);
Io = imopen(I,se);


Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);



Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);


fgm = imregionalmax(Iobrcbr);


I2 = imfuse(I, fgm);


se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);

fgm4 = bwareaopen(fgm3,20);
I3 = imfuse(I,fgm4);




img = imsharpen(I);
img = medfilt2(img);

Canny = edge(img, 'canny');



bw = bwdist(Iobrcbr);

D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;


%gmag2 = imimposemin(gmag, bgm | fgm4);
%L = watershed(gmag2);

%imshow(L,[])

%labels = imdilate(L == 0,ones(3,3) + 2*bgm + 3*fgm4); 
%I4 = imfuse(I,labels);
%imshow(I4,[])



