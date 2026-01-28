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
    (-VARWH 0 -VARB)
    (-VARDH 0 -VARB)
    ( VARDH 0 -VARB)
    ( VARWH 0 -VARB)

    (-VARWH VARH -VARB)
    (-VARDH VARH -VARB)
    ( VARDH VARH -VARB)
    ( VARWH VARH -VARB)

    (-VARWH 0 VARB)
    (-VARDH 0 VARB)
    ( VARDH 0 VARB)
    ( VARWH 0 VARB)

    (-VARWH VARH VARB)
    (-VARDH VARH VARB)
    ( VARDH VARH VARB)
    ( VARWH VARH VARB)
);

blocks
(
    hex (0 1 5 4 8 9 13 12) (VARNDW VARNY 1) simpleGrading (VARGXi VARGYi 1)
    hex (1 2 6 5 9 10 14 13) (VARND VARNY 1) simpleGrading (1 VARGYi 1)
    hex (2 3 7 6 10 11 15 14) (VARNDW VARNY 1) simpleGrading (VARGX VARGYi 1)
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
            (5 13 14 6)
        );
    }
    outlet
    {
        type wall;
        faces
        (
            (8 0 1 9)
            (9 1 2 10)
            (10 2 3 11)
        );
    }
    walls
    {
        type wall;
        faces
        (
            (12 4 0 8)
            (7 15 11 3)
        );
    }
    top
    {
        type wall;
        faces
        (
            (13 5 4 12)
            (15 7 6 14)
        );
    }
    empties
    {
        type empty;
        faces
        (
            (13 12 8 9)
            (14 13 9 10)
            (15 14 10 11)
            (5 1 0 4)
            (6 2 1 5)
            (7 3 2 6)
        );
    }
);

mergePatchPairs
(
);
