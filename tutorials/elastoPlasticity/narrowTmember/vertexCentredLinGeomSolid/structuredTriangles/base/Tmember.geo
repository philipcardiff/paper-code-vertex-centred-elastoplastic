// Average mesh spacing
dx = 0.02;

// Bounding box width
L = 0.15;

// Hole radius
r = 0.005;

// Width
a = 0.05;

// Out of plane depth
d = 0.015;

//Mesh parameters
//Radius of hole boundary layer
rb = 0.0075;

//Second boundary layer around hole
r2 = 0.05; 

//Cell multiplier
cn = 1;

// Points
Point(1) = {L-a, 0, 0, dx};
Point(2) = {L, 0, 0, dx};
Point(3) = {L, L, 0, dx};
Point(4) = {0, L, 0, dx};
Point(5) = {0, L-a, 0, dx};
Point(6) = {L-a-r, L-a, 0, dx/8};
Point(7) = {L-a, L-a-r, 0, dx/8};
Point(8) = {L-a, L-a, 0, dx};
Point(9) = {L-a-rb, L-a, 0, dx/8};
Point(10) = {L-a, L-a-rb, 0, dx/8};
Point(11) = {L-a-rb*Cos(45*Pi/180), L-a+rb*Sin(45*Pi/180), 0, dx/8};
Point(12) = {L-a-r*Cos(45*Pi/180), L-a+r*Sin(45*Pi/180), 0, dx/8};
Point(13) = {L-a+rb*Cos(45*Pi/180), L-a+rb*Sin(45*Pi/180), 0, dx/8};
Point(14) = {L-a+r*Cos(45*Pi/180), L-a+r*Sin(45*Pi/180), 0, dx/8};
Point(15) = {L-a+r*Cos(45*Pi/180), L-a-r*Sin(45*Pi/180), 0, dx/8};
Point(16) = {L-a+rb*Cos(45*Pi/180), L-a-rb*Sin(45*Pi/180), 0, dx/8};
Point(17) = {L-a-r2, L-a, 0, dx};
Point(18) = {L-a-r2, L, 0, dx};
Point(19) = {L-a, L-a-r2, 0, dx};
Point(20) = {L, L-a-r2, 0, dx};
Point(21) = {L-a+r, L-a, 0, dx/8};
Point(22) = {L-a+rb, L-a, 0, dx/8};
Point(23) = {L, L-a, 0, dx};
Point(24) = {L-a,L-a+r, 0, dx/8};
Point(25) = {L-a, L-a+rb, 0, dx/8};
Point(26) = {L-a, L, 0, dx};

//Point(11) = {rb*Cos(45*Pi/180), rb*Sin(45*Pi/180), 0, dx};
//Point(10) = {rb*Cos(45*Pi/180), L, 0, dx};
//Point(11) = {L, rb*Sin(45*Pi/180), 0, dx};
//Point(12) = {r*Cos(45*Pi/180), r*Sin(45*Pi/180), 0, dx};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 20};
Line(3) = {20, 23};
Line(4) = {3, 26};
Line(5) = {18, 4};
Line(6) = {4, 5};
Line(7) = {5, 17};
Line(8) = {17, 9};
Line(9) = {6, 9};
Line(10) = {10, 7};
Line(11) = {10, 19};
Line(12) = {19, 1};
Line(13) = {19, 20};
Line(14) = {17, 18};
Line(15) = {11, 18};
Line(16) = {12, 11};
Line(17) = {14, 13};
Line(18) = {13, 3};
Line(19) = {16, 15};
Line(20) = {16, 20};
Circle(21) = {6, 8, 12};
Circle(22) = {11, 8, 25};
Circle(23) = {13, 8, 22};
Circle(24) = {15, 8, 7};
Circle(25) = {9, 8, 11};
Circle(26) = {12, 8, 24};
Circle(27) = {14, 8, 21};
Circle(28) = {16, 8, 10};
Line(29) = {21, 22};
Line(30) = {22, 23};
Line(31) = {24, 25};
Line(32) = {26, 25};
Circle(33) = {24, 8, 14};
Circle(34) = {21, 8, 15};
Circle(35) = {25, 8, 13};
Circle(36) = {22, 8, 16};
Line(37) = {26, 18};
Line(38) = {23, 3};

// Surface
Curve Loop(1) = {7, 14, 5, 6};
Plane Surface(1) = {1};
Curve Loop(2) = {8, 25, 15, -14};
Plane Surface(2) = {2};
Curve Loop(3) = {15, -37, 32, -22};
Plane Surface(3) = {3};
Curve Loop(4) = {32, 35, 18, 4};
Plane Surface(4) = {4};
Curve Loop(5) = {18, -38, -30, -23};
Plane Surface(5) = {5};
Curve Loop(6) = {30, -3, -20, -36};
Plane Surface(6) = {6};
Curve Loop(7) = {20, -13, -11, -28};
Plane Surface(7) = {7};
Curve Loop(8) = {2, -13, 12, 1};
Plane Surface(8) = {8};
Curve Loop(9) = {25, -16, -21, 9};
Plane Surface(9) = {9};
Curve Loop(10) = {22, -31, -26, 16};
Plane Surface(10) = {10};
Curve Loop(11) = {35, -17, -33, 31};
Plane Surface(11) = {11};
Curve Loop(12) = {23, -29, -27, 17};
Plane Surface(12) = {12};
Curve Loop(13) = {36, 19, -34, 29};
Plane Surface(13) = {13};
Curve Loop(14) = {28, 10, -24, -19};
Plane Surface(14) = {14};

//Force mapped meshing (triangles)
Transfinite Curve {7, -5, 2, -12} = (10*cn)+1 Using Progression 1;
Transfinite Curve {6, 1} = (5*cn)+1 Using Progression 1;
//Transfinite Curve {34, 23, 27, 33, 35, 22, 26, 25, 21, 24, 28, 36, 14, 37, 4, 38, 3, 13} = (5*cn)+1 Using Progression 1;
Transfinite Curve {11, 15, 8, 18, 32, 30, 20} = (12*cn)+1 Using Progression 1;
Transfinite Curve {3, 23, 22, 4, 26, 27} = (14*cn)+1 Using Progression 1;
Transfinite Curve {9, 16, 17, 19, 10, 29, 31} = cn+1 Using Progression 1;
Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

//Recombine surface to get quadrilaterals
//Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

// Create volume by extrusion
Extrude {0, 0, d} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5};    
  Surface{6}; Surface{7}; Surface{8}; Surface{9}; Surface{10};
  Surface{11}; Surface{12}; Surface{13}; Surface{14};
  Layers{cn*2};
  //Recombine;
}

Physical Surface("loading") = {213};
Physical Surface("front") = {60, 82, 104, 126, 170, 192, 214, 346, 324, 302, 280, 258, 236, 148};
Physical Surface("fixed") = {59};
Physical Surface("symmx") = {201, 161, 139};
Physical Surface("symmz") = {6, 5, 4, 3, 2, 1, 10, 9, 8, 7, 11, 12, 13, 14};
Physical Surface("tractionFree") = {47, 69, 345, 301, 267, 245, 223, 187, 209, 125, 95, 55, 235, 231, 253, 275, 297, 319, 341, 337};
Physical Volume("internal") = {6, 5, 4, 3, 2, 1, 7, 8, 9, 10, 11, 12, 13, 14};


Transfinite Curve {34, 23, 27, 33, 35, 22, 26, 25, 21, 24, 28, 36, 14, 37, 4, 38, 3, 13} = (5*cn)+1 Using Progression 1;

