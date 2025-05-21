/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2019 Philip Morris International

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "mixedAbsorbingFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::mixedAbsorbingFvPatchScalarField::Peclet() const
{
    const aerosolModel& aerosol =
        db().lookupObject<aerosolModel>("aerosolProperties");

    const label i = patch().index();

    vectorField W(aerosol.U().boundaryField()[i]);

    // In some cases the drift velocity may not yet be in the database. In that
    // situation, use U as W.

    if (db().foundObject<volVectorField>(driftVelocityName_))
    {
        const vectorField& V =
            db().lookupObject<volVectorField>(driftVelocityName_)
           .boundaryField()[i];

        W += V;
    }

    const scalarField& D =
        db().lookupObject<volScalarField>(diffusivityName_)
       .boundaryField()[i];

    // Peclet

    return mag(patch().delta())*mag(W)/max(D,SMALL);
}

Foam::tmp<Foam::scalarField>
Foam::mixedAbsorbingFvPatchScalarField::inwardFluxMask() const
{
    const aerosolModel& aerosol =
        db().lookupObject<aerosolModel>("aerosolProperties");

    const label i = patch().index();

    const scalarField& rho =
        aerosol.rho().boundaryField()[i];

    // In some cases the drift velocity may not yet be in the database. In that
    // situation, set V as zero.

    vectorField V(patch().size(), Zero);

    if (db().foundObject<volVectorField>(driftVelocityName_))
    {
        V =
            db().lookupObject<volVectorField>(driftVelocityName_)
           .boundaryField()[i];
    }

    const vectorField& Sf = aerosol.mesh().Sf().boundaryField()[i];

    const scalarField phiInertial(rho * (V & Sf));

    const scalarField phi
    (
        aerosol.phi().boundaryField()[i]
      - aerosol.drift().phiCorr().boundaryField()[i]
      + phiInertial
    );

    return 1.0 - pos0(phi);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedAbsorbingFvPatchScalarField::
mixedAbsorbingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    driftVelocityName_("Ur"),
    diffusivityName_("DDisp"),
    gamma_(1.0),
    clipInwardFlux_(true)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = Zero;
}

Foam::mixedAbsorbingFvPatchScalarField::
mixedAbsorbingFvPatchScalarField
(
    const mixedAbsorbingFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    driftVelocityName_(ptf.driftVelocityName_),
    diffusivityName_(ptf.diffusivityName_),
    gamma_(ptf.gamma_),
    clipInwardFlux_(ptf.clipInwardFlux_)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
}

Foam::mixedAbsorbingFvPatchScalarField::
mixedAbsorbingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    driftVelocityName_(dict.lookupOrDefault<word>("driftVelocity", "Ur")),
    diffusivityName_(dict.lookupOrDefault<word>("diffusivity", "DDisp")),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.0)),
    clipInwardFlux_(dict.lookupOrDefault<bool>("clipInwardFlux", true))
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = Zero;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
       fvPatchField<scalar>::operator=(0.0);
    }
}

Foam::mixedAbsorbingFvPatchScalarField::
mixedAbsorbingFvPatchScalarField
(
    const mixedAbsorbingFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    driftVelocityName_(tppsf.driftVelocityName_),
    diffusivityName_(tppsf.diffusivityName_),
    gamma_(tppsf.gamma_),
    clipInwardFlux_(tppsf.clipInwardFlux_)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
}

Foam::mixedAbsorbingFvPatchScalarField::
mixedAbsorbingFvPatchScalarField
(
    const mixedAbsorbingFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    driftVelocityName_(tppsf.driftVelocityName_),
    diffusivityName_(tppsf.diffusivityName_),
    gamma_(tppsf.gamma_),
    clipInwardFlux_(tppsf.clipInwardFlux_)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedAbsorbingFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}

void Foam::mixedAbsorbingFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)

{
    mixedFvPatchScalarField::rmap(ptf, addr);
}

void Foam::mixedAbsorbingFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Peclet number

    const scalarField Pe(this->Peclet());

    // Set the value fraction

    this->valueFraction() = 1.0/(pow(Pe,gamma_) + 1.0);

    if (clipInwardFlux_)
    {
        // Clip inward flux by setting homogeneous Dirichlet

        const scalarField mask(this->inwardFluxMask());

        this->valueFraction() =
            this->valueFraction()*(1.0 - mask) + mask;
    }

    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::mixedAbsorbingFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntry("value", os);

    os.writeEntryIfDifferent<word>("driftVelocity", "Ur", driftVelocityName_);
    os.writeEntryIfDifferent<word>("diffusivity", "DDisp", diffusivityName_);
    os.writeEntryIfDifferent<scalar>("gamma", 1.0, gamma_);
    os.writeEntryIfDifferent<bool>("clipInwardFlux", true, clipInwardFlux_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixedAbsorbingFvPatchScalarField
    );
}

// ************************************************************************* //
