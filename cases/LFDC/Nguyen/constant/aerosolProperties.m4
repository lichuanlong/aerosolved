FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      aerosolProperties;
}

aerosolModel    VARMODEL;

diameter
{
    min     1E-9;
    max     1E-4;
}

twoMomentLogNormalAnalyticalCoeffs
{
    sigma   1.2;
}

fixedSectionalCoeffs
{
    distribution
    {
        type    logarithmic;
        yMin    VARYMIN;
        yMax    VARYMAX;
        N       VARN;
    }

    interpolation
    {
        type    twoMoment;
    }

    rescale     true;

    initFromPatch inlet;
}

submodels
{
    condensation
    {
        type            coupled;

        KelvinEffect    true;

        heatOfVaporization
        {
            active          false;
        }
    }

    nucleation
    {
        type        coupled;
        tolerance   1E-7;
    }

    coalescence
    {
        type        none;
    }

    driftFluxModel
    {
        continuousDiffusion
        {
            type        none;
        }

        dispersedDiffusion
        {
            type        none;
        }

        dispersedInertialDrift
        {
            type        none;
        }
    }
}
