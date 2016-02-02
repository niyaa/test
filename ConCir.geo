// Gmsh project created on Tue Jan 26 12:16:46 2016
DefineConstant[ lc = { 4, Path "Gmsh/Parameters"}];
Point(1) = {0, 0, 0,  lc};
Point(2) = {1, 0, 0,  lc};
Point(3) = {-1, 0, 0,  lc};
Point(4) = {4, 0, 0,  lc};
Point(5) = {-4, 0, 0, lc};
Circle(1) = {5, 1, 4};
Circle(2) = {4, 1, 5};
Circle(3) = {3, 1, 2};
Circle(4) = {2, 1, 3};
Line Loop(5) = {2, 1};
Line Loop(6) = {4, 3};
Plane Surface(7) = {5, 6};
Physical Line("In") = {1, 2};
Physical Line("Cyl") = {4, 3};
Physical Surface("Domain") = {7};


