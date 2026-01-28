#include "zeroGradientDepositionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    zeroGradientMode_("V"),
    clipInwardFlux_(true)
{}


Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    zeroGradientMode_
    (
        dict.lookupOrDefault<word>("zeroGradientMode", "V")
    ),
    clipInwardFlux_
    (
        dict.lookupOrDefault<bool>("clipInwardFlux", true)
    )
{
    if (zeroGradientMode_ != "V" && zeroGradientMode_ != "Ud")
    {
        FatalErrorInFunction
            << "zeroGradientMode should be either V or Ud"
            << abort(FatalError);
    }
}


Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const zeroGradientDepositionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    zeroGradientMode_(ptf.zeroGradientMode_),
    clipInwardFlux_(ptf.clipInwardFlux_)
{}


Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const zeroGradientDepositionVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    zeroGradientMode_(ptf.zeroGradientMode_),
    clipInwardFlux_(ptf.clipInwardFlux_)
{}


Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const zeroGradientDepositionVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    zeroGradientMode_(ptf.zeroGradientMode_),
    clipInwardFlux_(ptf.clipInwardFlux_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroGradientDepositionVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField phip(patch().size(), 0.0);

    if (db().foundObject<aerosolModel>("aerosolProperties"))
    {
        const aerosolModel& aerosol =
            db().lookupObject<aerosolModel>("aerosolProperties");

        // Apply zero-gradient to either Ud or V

        vectorField V(patch().size());

        if (zeroGradientMode_ == "Ud")
        {
            // Apply zero-gradient to Ud and then reconstruct V from that. The
            // relation assumes that the corrective velocity at the patch is
            // equal to the cell center one

            V =
                aerosol.U().boundaryField()[patch().index()].patchInternalField()
              - aerosol.U().boundaryField()[patch().index()]
              + this->patchInternalField();
        }
        else
        {
            // Otherwise, initialize V as zero-gradient directly

            V = this->patchInternalField();
        }

        if (clipInwardFlux_)
	    {
       	    // Compute the mass flux associated with V

            const scalarField phiV
        	(
        	    (V & patch().nf())*patch().magSf()
        	  * aerosol.rho().boundaryField()[patch().index()]
        	);

        	// Compute the dispersed phase flux from Ud = U + V - VCorr

        	const scalarField phiUd
        	(
        	    aerosol.phi().boundaryField()[patch().index()]
        	  - aerosol.drift().phiCorr().boundaryField()[patch().index()]
        	  + phiV
        	);

        	// For faces which have a flux into the fluid domain, set a relative
        	// velocity such that the dispersed velocity is exactly zero

        	const scalarField phiZero
        	(
        	    aerosol.drift().phiCorr().boundaryField()[patch().index()]
        	  - aerosol.phi().boundaryField()[patch().index()]
        	);

            operator==
            (
                V*pos0(phiUd)
              + (1.0-pos0(phiUd))*phiZero*patch().nf()
              / (patch().magSf()*aerosol.rho().boundaryField()[patch().index()])
            );
	    }
        else
	    {
	        operator==(V);
	    }
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::zeroGradientDepositionVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

    os.writeEntryIfDifferent<word>("zeroGradientMode", "V", zeroGradientMode_);
    os.writeEntryIfDifferent<bool>("clipInwardFlux", true, clipInwardFlux_);

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       zeroGradientDepositionVelocityFvPatchVectorField
   );
}

// ************************************************************************* //
