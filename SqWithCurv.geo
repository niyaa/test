// Gmsh project created on Tue Jan 26 22:45:40 2016
DefineConstant[ lc = { 1, Path "Gmsh/Parameters"}];
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {1, 1, 0, lc};
Point(5) = {0, 1, 0, lc};
Point(6) = {1, 0.5, 0, lc};
Line(1) = {5, 1};
Line(2) = {1, 2};
Line(3) = {5, 3};
Line(4) = {2, 6};
Delete {
  Line{4};
}
Circle(4) = {2, 6, 3};
Line Loop(5) = {1, 2, 4, -3};
Plane Surface(6) = {5};
Physical Line("In") = {1};
Physical Line("Wall") = {3, 2};
Physical Line("Cyl") = {4};
Physical Surface("Domain") = {6};
