FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Butanol.continuous;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    inlet
    {
        type        saturatedMixture;
        value       $internalField;
        T           VARTSAT;
    }

    outlet
    {
        type        inletOutlet;
        inletValue  $internalField;
        value       $internalField;
    }

    wallPreheater
    {
        type        saturatedMixture;
        value       $internalField;
    }

    wallCondenser
    {
        type        saturatedMixture;
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
