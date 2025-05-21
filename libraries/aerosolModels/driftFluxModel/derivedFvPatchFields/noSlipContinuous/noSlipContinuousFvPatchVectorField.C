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

#include "noSlipContinuousFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noSlipContinuousFvPatchVectorField::
noSlipContinuousFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::noSlipContinuousFvPatchVectorField::
noSlipContinuousFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::noSlipContinuousFvPatchVectorField::
noSlipContinuousFvPatchVectorField
(
    const noSlipContinuousFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::noSlipContinuousFvPatchVectorField::
noSlipContinuousFvPatchVectorField
(
    const noSlipContinuousFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf)
{}


Foam::noSlipContinuousFvPatchVectorField::
noSlipContinuousFvPatchVectorField
(
    const noSlipContinuousFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::noSlipContinuousFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (db().foundObject<aerosolModel>("aerosolProperties"))
    {
        const aerosolModel& aerosol =
            db().lookupObject<aerosolModel>("aerosolProperties");

        const aerosolThermo& thermo = aerosol.thermo();

        const speciesTable& dispSpecies = thermo.dispSpecies();

        // Collect all dispersed phase effective fluxes

        scalarField phiEff(patch().size(), 0.0);

        forAll(dispSpecies, j)
        {
            phiEff +=
                patch().lookupPatchField<surfaceScalarField, vector>
                (
                    IOobject::groupName
                    (
                        "phiEff",
                        thermo.Z()[j].name()
                    )
                );
        }

        // Set the velocity to be equal to the dispersed phase velocity, which
        // implies that the continuous phase velocity is zero

        const scalarField rho(aerosol.thermo().rho(patch().index()));

        vectorField::operator=
        (
            phiEff/(this->patch().magSf()*rho)
          * patch().nf()
        );
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::noSlipContinuousFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        noSlipContinuousFvPatchVectorField
    );
}

// ************************************************************************* //
