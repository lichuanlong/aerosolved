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

#include "twoMomentSubGridDepositionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "aerosolModel.H"
#include "subGridDepositionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoMomentSubGridDepositionVelocityFvPatchVectorField::
twoMomentSubGridDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    subGridDepositionVelocityFvPatchVectorField(p, iF)
{}


Foam::twoMomentSubGridDepositionVelocityFvPatchVectorField::
twoMomentSubGridDepositionVelocityFvPatchVectorField
(
    const twoMomentSubGridDepositionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    subGridDepositionVelocityFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::twoMomentSubGridDepositionVelocityFvPatchVectorField::
twoMomentSubGridDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    subGridDepositionVelocityFvPatchVectorField(p, iF, dict)
{}


Foam::twoMomentSubGridDepositionVelocityFvPatchVectorField::
twoMomentSubGridDepositionVelocityFvPatchVectorField
(
    const twoMomentSubGridDepositionVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    subGridDepositionVelocityFvPatchVectorField(fcvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::twoMomentSubGridDepositionVelocityFvPatchVectorField::d() const
{
    const aerosolModel& aerosol =
        db().lookupObject<aerosolModel>("aerosolProperties");

    return aerosol.meanDiameter(5.0,3.0,patch().index());
}

void Foam::twoMomentSubGridDepositionVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    subGridDepositionVelocityFvPatchVectorField::evaluate(commsType);
}

void Foam::twoMomentSubGridDepositionVelocityFvPatchVectorField::write
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
        twoMomentSubGridDepositionVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
