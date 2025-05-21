FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      T;
}

dimensions      [0 0 0 1 0 0 0];

internalField   uniform VARTPREH;

boundaryField
{
    inlet
    {
        type        fixedValue;
        value       uniform VARTPREH;
    }

    outlet
    {
        type        inletOutlet;
        inletValue  $internalField;
        value       $internalField;
    }

    wallPreheater
    {
        type        fixedValue;
        value       uniform VARTPREH;
    }

    wallCondenser
    {
        type        fixedValue;
        value       uniform VARTCOND;
    }

    front
    {
        type        wedge;
    }

    back
    {
        type        wedge;
    }

    axis
    {
        type        empty;
    }
}
