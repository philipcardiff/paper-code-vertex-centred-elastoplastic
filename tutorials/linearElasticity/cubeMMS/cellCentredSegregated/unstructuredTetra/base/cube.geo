// Average mesh spacing
dx = 0.5;

// Bounding box width
L = 1; 

//Cell multiplier
cn = 1;

// Points
Point(1) = {0, 0, 0, dx};
Point(2) = {L, 0, 0, dx};
Point(3) = {0, L, 0, dx};
Point(4) = {L, L, 0, dx};

// Lines
Line(1) = {3, 1};
Line(2) = {1, 2};
Line(3) = {2, 4};
Line(4) = {4, 3};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//Force mapped meshing (triangles)
//Transfinite Curve {2, 3, 4, 1} = cn+1 Using Progression 1;
//Transfinite Surface {1};

//Recombine surface to get quadrilaterals
// Recombine Surface {1};

// Create volume by extrusion
Extrude {0, 0, L} {
  Surface{1};
//  Layers{cn};
//  Recombine;
}


// Physical groups
Physical Surface("right") = {21};
Physical Surface("left") = {13};
Physical Surface("back") = {25};
Physical Surface("front") = {17};
Physical Surface("top") = {26};
Physical Surface("bottom") = {1};
Physical Volume("internal") = {1};
