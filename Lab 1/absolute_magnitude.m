function mag = absolute_magnitude(Gx,Gy)
% Compute the approximate magnitude of the gradients
% Input: Gx, Gy - the gradients in the x and y directions
% Output: mag - the approximate magnitude
mag = abs(Gx) + abs(Gy);
show_image(mag);
show_image(mag>20);
show_image(mag>40);
show_image(mag>60);
show_image(mag>80);