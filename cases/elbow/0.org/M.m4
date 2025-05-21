FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      M;
}

dimensions      [-1 0 0 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    inlet
    {
        type        VARMINLETBC;
        CMD         1E-6;
        sigma       0.5;
        value       $internalField;
    }

    outlet
    {
        type        inletOutlet;
        inletValue  $internalField;
        value       $internalField;
    }

    walls
    {
        type        VARMWALLBC;
        value       $internalField;
    }

    empties
    {
        type        empty;
    }
}
