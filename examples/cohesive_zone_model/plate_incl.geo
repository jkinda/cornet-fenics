d = 0.015;
R = 0.25;
Point(1) = {0., 0., 0., d};
Point(2) = {0.4-R, 0., 0., d};
Point(3) = {0.4, 0., 0., d};
Point(4) = {0.4+R, 0., 0., d};
Point(5) = {1., 0., 0., d};
Point(6) = {1., 0.5, 0., d};
Point(7) = {0.6+R, 0.5, 0., d};
Point(8) = {0.6, 0.5, 0., d};
Point(9) = {0.6-R, 0.5, 0., d};
Point(10) = {0., 0.5, 0., d};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 1};
Circle(11) = {4,3,2};
Circle(12) = {9,8,7};
//+
Line Loop(1) = {4, 5, 6, -12, 9, 10, 1, -11};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {11, 2, 3};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {12, 7, 8};
//+
Plane Surface(3) = {3};

Physical Line(1) = {10};
Physical Line(2) = {5};
Physical Line(3) = {11, 12};

Physical Surface(1) = {1};
Physical Surface(2) = {2,3};
