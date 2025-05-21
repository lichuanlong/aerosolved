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
    (-VARWH 0 -VARLH)
    ( VARWH 0 -VARLH)
    ( VARWH VARH -VARLH)
    (-VARWH VARH -VARLH)

    (-VARWH 0 VARLH)
    ( VARWH 0 VARLH)
    ( VARWH VARH VARLH)
    (-VARWH VARH VARLH)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (VARNW VARNH VARNL)

    simpleGrading
    (
        (
            (0.4 0.35 0.5)
            (0.2 0.3 1)
            (0.4 0.35 2)
        )
        1
        (
            (0.4 0.35 0.5)
            (0.2 0.3 1)
            (0.4 0.35 2)
        )
    )
);

edges
(
);

boundary
(
    top
    {
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }
    outlet
    {
        type wall;
        faces
        (
            (4 0 1 5)
        );
    }
    walls
    {
        type wall;
        faces
        (
            (7 3 0 4)
            (2 6 5 1)
            (7 4 5 6)
            (3 2 1 0)
        );
    }
);

mergePatchPairs
(
);
