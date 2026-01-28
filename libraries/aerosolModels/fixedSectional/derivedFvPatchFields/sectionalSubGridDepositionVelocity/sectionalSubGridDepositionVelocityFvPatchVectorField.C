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

#include "sectionalSubGridDepositionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "aerosolModel.H"
#include "subGridDepositionModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label
Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::sectionIndex()
const
{
    return readLabel(IStringStream(internalField().group())()) - 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::
sectionalSubGridDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    subGridDepositionVelocityFvPatchVectorField(p, iF)
{}


Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::
sectionalSubGridDepositionVelocityFvPatchVectorField
(
    const sectionalSubGridDepositionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    subGridDepositionVelocityFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::
sectionalSubGridDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    subGridDepositionVelocityFvPatchVectorField(p, iF, dict)
{}

Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::
sectionalSubGridDepositionVelocityFvPatchVectorField
(
    const sectionalSubGridDepositionVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    subGridDepositionVelocityFvPatchVectorField(fcvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::d() const
{
    const fixedSectionalSystem& system =
        db().lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

    const aerosolModel& aerosol = system.aerosol();
    const aerosolThermo& thermo = aerosol.thermo();

    const scalarField rhod(thermo.thermoDisp().rho(patch().index()));

    const section& sec = system.distribution()[this->sectionIndex()];

    return sec.d(rhod);
}

void Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    subGridDepositionVelocityFvPatchVectorField::evaluate(commsType);
}

void Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    subGridDepositionVelocityFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        sectionalSubGridDepositionVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
