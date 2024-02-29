shakey = read_image('','shakey.150.gif');
%{show_image(shakey)
load sobel 
shakey_sobelX = conv2(shakey,sobelX,'valid');
shakey_sobelY = conv2(shakey,sobelY,'valid');
%{show_image(abs(shakey_sobelX)>5)
magnitude_sobel = magnitude(shakey_sobelX,shakey_sobelY);
load roberts
shakey_robertsA = conv2(shakey,robertsA,'valid');
shakey_robertsB = conv2(shakey,robertsB,'valid');
%{magnitude_roberts = magnitude(shakey_robertsA, shakey_robertsB);
absolute_mag = absolute_magnitude(shakey_sobelX,shakey_sobelY);