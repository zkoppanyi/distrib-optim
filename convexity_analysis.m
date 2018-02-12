clear variables;

syms x1 x2 y1 y2 r
f = (x1 - x2)^2 + (y1 - y2)^2 - r^2;

h(1,1) = diff(diff(f, x1), x1);
h(1,2) = diff(diff(f, x1), y1);
h(1,3) = diff(diff(f, x1), x2);
h(1,4) = diff(diff(f, x1), y2);

h(2,1) = diff(diff(f, y1), x1);
h(2,2) = diff(diff(f, y1), y1);
h(2,3) = diff(diff(f, y1), x2);
h(2,4) = diff(diff(f, y1), y2);

h(3,1) = diff(diff(f, x2), x1);
h(3,2) = diff(diff(f, x2), y1);
h(3,3) = diff(diff(f, x2), x2);
h(3,4) = diff(diff(f, x2), y2);

h(4,1) = diff(diff(f, y2), x1);
h(4,2) = diff(diff(f, y2), y1);
h(4,3) = diff(diff(f, y2), x2);
h(4,4) = diff(diff(f, y2), y2);

H = [2 0 -2 0;0,  2,  0, -2; -2,  0,  2,  0;0, -2,  0,  2];

%%
syms x1 x2 x3 x4 y1 y2 y3 y4 r
f = ((x1 - x2)^2 + (y1 - y2)^2 - r^2)^2 + ((x1 - x3)^2 + (y1 - y3)^2 - r^2)^2 + ((x1 - x4)^2 + (y1 - y4)^2 - r^2)^2;

h2(1,1) = diff(diff(f, x1), x1);
h2(1,2) = diff(diff(f, x1), y1);

h2(2,1) = diff(diff(f, y1), x1);
h2(2,2) = diff(diff(f, y1), y1);

H = [ 6,  0, -2,  0;
      0,  6,  0, -2;
     -2,  0,  2,  0;
      0, -2,  0,  2];
  
  
[ 4*(x1 - x2)^2 + 4*(x1 - x3)^2 + 4*(x1 - x4)^2 + 4*(y1 - y2)^2 + 4*(y1 - y3)^2 + 4*(y1 - y4)^2 + 2*(2*x1 - 2*x2)^2 + 2*(2*x1 - 2*x3)^2 + 2*(2*x1 - 2*x4)^2 - 12*r^2,                                                                      2*(2*x1 - 2*x2)*(2*y1 - 2*y2) + 2*(2*x1 - 2*x3)*(2*y1 - 2*y3) + 2*(2*x1 - 2*x4)*(2*y1 - 2*y4)]
[                                                                      2*(2*x1 - 2*x2)*(2*y1 - 2*y2) + 2*(2*x1 - 2*x3)*(2*y1 - 2*y3) + 2*(2*x1 - 2*x4)*(2*y1 - 2*y4), 4*(x1 - x2)^2 + 4*(x1 - x3)^2 + 4*(x1 - x4)^2 + 4*(y1 - y2)^2 + 4*(y1 - y3)^2 + 4*(y1 - y4)^2 + 2*(2*y1 - 2*y2)^2 + 2*(2*y1 - 2*y3)^2 + 2*(2*y1 - 2*y4)^2 - 12*r^2]
  