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
        name    preheater;
        type    faceSet;
        action  new;
        source  boxToFace;

        sourceInfo
        {
            min     (0 -99 VARZ2);
            max     (VARLPREH 99 99);
        }
    }
);
