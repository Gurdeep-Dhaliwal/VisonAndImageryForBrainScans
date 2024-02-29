function m = magnitude(x,y)
%{m = hypot(x,y);}%
m = sqrt(x.^2+y.^2);
show_image(m);
show_image(m>20);
show_image(m>40);
show_image(m>60);
show_image(m>80);
