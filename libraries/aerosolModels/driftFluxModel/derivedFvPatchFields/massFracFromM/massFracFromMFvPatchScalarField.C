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

#include "massFracFromMFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

#include "fixedSectionalSystem.H"
#include "aerosolModel.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massFracFromMFvPatchScalarField::
massFracFromMFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    absorbingMFvPatchScalarField(p, iF)
{}


Foam::massFracFromMFvPatchScalarField::
massFracFromMFvPatchScalarField
(
    const massFracFromMFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    absorbingMFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::massFracFromMFvPatchScalarField::
massFracFromMFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    absorbingMFvPatchScalarField(p, iF, dict)
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
    else
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
}


Foam::massFracFromMFvPatchScalarField::
massFracFromMFvPatchScalarField
(
    const massFracFromMFvPatchScalarField& tppsf
)
:
    absorbingMFvPatchScalarField(tppsf)
{}


Foam::massFracFromMFvPatchScalarField::
massFracFromMFvPatchScalarField
(
    const massFracFromMFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    absorbingMFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::massFracFromMFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (db().foundObject<fixedSectionalSystem>("fixedSectionalSystem"))
    {
        const fixedSectionalSystem& system =
            db().lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

        const aerosolModel& aerosol = system.aerosol();
        const aerosolThermo& thermo = aerosol.thermo();
        const speciesTable& activeSpecies = thermo.activeSpecies();

        scalarField alpha(patch().size(), 0.0);
        scalarField alphaf(patch().size(), 0.0);

        forAll(activeSpecies, j)
        {
            alpha +=
                thermo.Z()[j].boundaryField()[patch().index()]
               .patchInternalField();
        }

        forAll(system.distribution(), i)
        {
            alphaf +=
                system.distribution()[i].x()
              * system.distribution()[i].M()
               .boundaryField()[patch().index()];
        }

        const scalarField frac(this->patchInternalField()/max(alpha,SMALL));

        operator==(frac*alphaf);

        fixedValueFvPatchScalarField::updateCoeffs();
    }
    else
    {
        absorbingMFvPatchScalarField::updateCoeffs();
    }
}


void Foam::massFracFromMFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        massFracFromMFvPatchScalarField
    );
}

// ************************************************************************* //
