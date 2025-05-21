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

#include "noSlipInertFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noSlipInertFvPatchVectorField::
noSlipInertFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::noSlipInertFvPatchVectorField::
noSlipInertFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::noSlipInertFvPatchVectorField::
noSlipInertFvPatchVectorField
(
    const noSlipInertFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::noSlipInertFvPatchVectorField::
noSlipInertFvPatchVectorField
(
    const noSlipInertFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf)
{}


Foam::noSlipInertFvPatchVectorField::
noSlipInertFvPatchVectorField
(
    const noSlipInertFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::noSlipInertFvPatchVectorField::updateCoeffs()
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

        const speciesTable& contSpecies = thermo.contSpecies();
        const speciesTable& dispSpecies = thermo.dispSpecies();

        // Collect all non-inert effective fluxes

        scalarField phiEff(patch().size(), 0.0);

        forAll(contSpecies, j)
        {
            if (contSpecies[j] != thermo.inertSpecie())
            {
                phiEff +=
                    patch().lookupPatchField<surfaceScalarField, vector>
                    (
                        IOobject::groupName
                        (
                            "phiEff",
                            thermo.Y()[j].name()
                        )
                    );
            }
        }

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

        // Set the velocity to be equal to the sum of non-inert velocities,
        // which implies that the inert velocity is zero

        const scalarField rho(aerosol.thermo().rho(patch().index()));

        vectorField::operator=
        (
            phiEff/(this->patch().magSf()*rho)
          * patch().nf()
        );
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::noSlipInertFvPatchVectorField::write(Ostream& os) const
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
        noSlipInertFvPatchVectorField
    );
}

// ************************************************************************* //
