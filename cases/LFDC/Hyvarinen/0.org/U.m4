FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (VARU 0 0);

boundaryField
{
    inlet
    {
        type        fixedValue;
        value       uniform (VARU 0 0);
    }

    outlet
    {
        type        zeroGradient;
    }

    wallPreheater
    {
        type        noSlip;
    }

    wallCondenser
    {
        type        noSlip;
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
