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

#include "massFracFromSectionalFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

#include "fixedSectionalSystem.H"
#include "aerosolModel.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massFracFromSectionalFvPatchScalarField::
massFracFromSectionalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::massFracFromSectionalFvPatchScalarField::
massFracFromSectionalFvPatchScalarField
(
    const massFracFromSectionalFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::massFracFromSectionalFvPatchScalarField::
massFracFromSectionalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false)
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


Foam::massFracFromSectionalFvPatchScalarField::
massFracFromSectionalFvPatchScalarField
(
    const massFracFromSectionalFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::massFracFromSectionalFvPatchScalarField::
massFracFromSectionalFvPatchScalarField
(
    const massFracFromSectionalFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::massFracFromSectionalFvPatchScalarField::updateCoeffs()
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
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::massFracFromSectionalFvPatchScalarField::write(Ostream& os) const
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
        massFracFromSectionalFvPatchScalarField
    );
}

// ************************************************************************* //
