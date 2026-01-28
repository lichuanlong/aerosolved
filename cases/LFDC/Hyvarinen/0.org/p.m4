FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      p;
}

dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform VARP;

boundaryField
{
    inlet
    {
        type        fixedFluxPressure;
        value       $internalField;
    }

    outlet
    {
        type        totalPressure;
        p0          $internalField;
    }

    wallPreheater
    {
        type        fixedFluxPressure;
        value       $internalField;
    }

    wallCondenser
    {
        type        fixedFluxPressure;
        value       $internalField;
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
