FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

/*

Z6    o-----o
      |     |
      |     |
      |     |
      |     |
      |     |
      |     |
      | VIII|
      |     |
      |     |
      |     |
      |     |
      |     |
      |     |
Z5    o-----o
      |      \
      |       \
Z4    |        o
      |  VII  / \
      |      /   \
Z3    |     /     o--------o
      |    /      |        |
      |   /       |        |
      |  /   V    |   VI   |
      | /         |        |
      |/          |        |
Z2    o-----------o--------o
      |           |        |
      |     III   |   IV   |
      |           |        |
Z1    o-----------o--------o
      |     I     |   II   |
Z0    o-----------o--------o
     X0    X1 X2 X3       X4

*/

scale 0.001;

vertices
(
    (VARX0 VARY0 VARZ0)
    (VARX0 VARY0 VARZ1)
    (VARX0 VARY0 VARZ2)
    (VARX0 VARY0 VARZ5)
    (VARX0 VARY0 VARZ6)

    (VARX3 VARY3 VARZ0)
    (VARX3 VARY3 VARZ1)
    (VARX3 VARY3 VARZ2)
    (VARX3 VARY3 VARZ3)
    (VARX2 VARY2 VARZ4)
    (VARX1 VARY1 VARZ5)
    (VARX1 VARY1 VARZ6)

    (VARX4 VARY4 VARZ0)
    (VARX4 VARY4 VARZ1)
    (VARX4 VARY4 VARZ2)
    (VARX4 VARY4 VARZ3)

    (VARX3 -VARY3 VARZ0)
    (VARX3 -VARY3 VARZ1)
    (VARX3 -VARY3 VARZ2)
    (VARX3 -VARY3 VARZ3)
    (VARX2 -VARY2 VARZ4)
    (VARX1 -VARY1 VARZ5)
    (VARX1 -VARY1 VARZ6)

    (VARX4 -VARY4 VARZ0)
    (VARX4 -VARY4 VARZ1)
    (VARX4 -VARY4 VARZ2)
    (VARX4 -VARY4 VARZ3)
);

blocks
(
    hex (0 16 5 0 1 17 6 1) (VARNX1 1 VARNBL) simpleGrading (1 1 VARGBL)
    hex (16 23 12 5 17 24 13 6) (VARNX2 1 VARNBL) simpleGrading (1 1 VARGBL)
    hex (1 17 6 1 2 18 7 2) (VARNX1 1 VARNZ1) simpleGrading (1 1 1)
    hex (17 24 13 6 18 25 14 7) (VARNX2 1 VARNZ1) simpleGrading (1 1 1)
    hex (2 18 7 2 20 19 8 9) (VARNX1 1 VARNZ2) simpleGrading (1 1 1)
    hex (18 25 14 7 19 26 15 8) (VARNX2 1 VARNZ2) edgeGrading (1 1 VARGOUTLET VARGOUTLET 1 1 1 1 1 1 1 1)
    hex (2 20 9 2 3 21 10 3) (VARNZ2 1 VARNX1) simpleGrading (1 1 1)
    hex (3 21 10 3 4 22 11 4) (VARNZ2 1 VARNPIPE) edgeGrading (1 1 1 1 1 1 1 1 VARGPIPE VARGPIPE2 VARGPIPE2 VARGPIPE)
);

edges
(
    arc 9 10 (VARXC0 VARYC0 VARZC0)
    arc 20 21 (VARXC0 -VARYC0 VARZC0)

    arc 8 9 (VARXC1 VARYC1 VARZC1)
    arc 19 20 (VARXC1 -VARYC1 VARZC1)
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (4 22 11 4)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (19 26 15 8)
        );
    }

    lower
    {
        type wall;
        faces
        (
             (0 5 16 0)
             (5 12 23 16)
        );
    }

    walls
    {
        type wall;
        faces
        (
            (23 12 13 24)
            (24 13 14 25)
            (25 14 15 26)
            (19 8 9 20)
            (20 9 10 21)
            (21 10 11 22)
        );
    }

    rear
    {
        type wedge;
        faces
        (
            (5 0 1 6)
            (12 5 6 13)
            (6 1 2 7)
            (13 6 7 14)
            (7 2 9 8)
            (14 7 8 15)
            (9 2 3 10)
            (10 3 4 11)
         );
    }

    front
    {
        type wedge;
        faces
        (
            (0 16 17 1)
            (16 23 24 17)
            (1 17 18 2)
            (17 24 25 18)
            (2 18 19 20)
            (18 25 26 19)
            (2 20 21 3)
            (3 21 22 4)
         );
    }

    axis
    {
        type empty;
        faces
        (
            (0 1 1 0)
            (1 2 2 1)
            (2 3 3 2)
            (3 4 4 3)
        );
    }
);