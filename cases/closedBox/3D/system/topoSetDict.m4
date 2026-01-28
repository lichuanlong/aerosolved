FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}

actions
(
    {
        name    inlet;
        type    faceSet;
        action  new;
        source  cylinderToFace;
        p1      (0 VARP1 0);
        p2      (0 VARP2 0);
        radius  VARR;

        // For the foundation OpenFOAM version:

        sourceInfo
        {
            p1      (0 VARP1 0);
            p2      (0 VARP2 0);
            radius  VARR;
        }
    }
);
