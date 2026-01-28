FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      setFieldsDict;
}

regions
(
    boxToCell
    {
        box (0 -1 0) (VARLPREH 1 1);

        fieldValues
        (
            volScalarFieldValue T VARTPREH
        );
    }

    boxToCell
    {
        box (VARLPREH -1 0) (VARL 1 1);

        fieldValues
        (
            volScalarFieldValue T VARTCOND
        );
    }
);
