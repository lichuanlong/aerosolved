FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

convertToMeters 1;

vertices
(
    (-VARRX  0 -VARRY)
    (-VARRIX 0 -VARRIY)
    ( 0      0  0)

    ( 0      VARH  0)
    (-VARRIX VARH -VARRIY)
    (-VARRX  VARH -VARRY)

    (-VARRX  0 VARRY)
    (-VARRIX 0 VARRIY)

    (-VARRIX VARH VARRIY)
    (-VARRX  VARH VARRY)


);

blocks
(
    hex (0 1 4 5 6 7 8 9) (VARNRRI VARNY 1) simpleGrading (VARGRi VARGYi 1)
    hex (1 2 3 4 7 2 3 8) (VARNRI VARNY 1) simpleGrading (1 VARGYi 1)
);

edges
(
);

boundary
(
    inlet
    {
        type wall;
        faces
        (
            (3 4 8 3)
        );
    }
    outlet
    {
        type wall;
        faces
        (
            (7 6 0 1)
            (2 7 1 2)
        );
    }
    walls
    {
        type wall;
        faces
        (
            (9 5 0 6)
        );
    }
    top
    {
        type wall;
        faces
        (
            (4 5 9 8)
        );
    }
    front
    {
        type wedge;
        faces
        (
            (8 9 6 7)
            (3 8 7 2)
        );
    }
    back
    {
        type wedge;
        faces
        (
            (5 4 1 0)
            (4 3 2 1)
        );
    }
    axis
    {
        type symmetry;
        faces
        (
            (3 2 2 3)
        );
    }
);

mergePatchPairs
(
);
