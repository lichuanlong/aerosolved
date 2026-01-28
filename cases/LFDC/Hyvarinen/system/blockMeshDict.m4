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
    (0      0       0)
    (VARL   0       0)
    (0     -VARY    VARZ)
    (VARL  -VARY    VARZ)
    (VARL   VARY    VARZ)
    (0      VARY    VARZ)
);

blocks
(
    hex (0 1 1 0 2 3 4 5) (VARNX 1 VARNZ) simpleGrading (1 1 1)
);


edges
(
);

boundary
(
    axis
    {
        type empty;
        faces
        (
            (0 1 1 0)
        );
    }

    inlet
    {
        type patch;
        faces
        (
            (0 2 5 0)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (1 3 4 1)
        );
    }

    wallCondenser
    {
        type wall;
        faces
        (
            (2 3 4 5)
        );
    }

    back
    {
        type wedge;
        faces
        (
            (0 1 4 5)
        );
    }

    front
    {
        type wedge;
        faces
        (
            (0 1 3 2)
        );
    }
);

mergePatchPairs
(
);
