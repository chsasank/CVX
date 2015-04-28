l = 1;
delta = 0.01;
[X,Y,W]=meshgrid(0:delta:l,0:delta:l,-l*2:delta:l*2);
% isosurface(X,Y,W,(X.*Y-W.^2).*(X-l).*(Y-l),0)



V = ones(size(X));
for theta = -pi:0.05:pi
    s = [cos(theta) ;sin(theta)];
    V = V.*( X*cos(theta)^2 + 2*W*cos(theta)*sin(theta) + Y*sin(theta)^2 >=0);
end

figure
isosurface(X,Y,W,V,0)
xlabel('x');ylabel('y');zlabel('w')
title('isosurface of the set from intersection of halfplanes')

figure
isosurface(X,Y,W,(X.*Y-W.^2),0)
xlabel('x');ylabel('y');zlabel('w')
title('isosurface plot of the set from determinant >= 0')