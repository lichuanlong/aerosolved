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

#include "sectionalMixedAbsorbingFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedSectional.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label
Foam::sectionalMixedAbsorbingFvPatchScalarField::sectionIndex() const
{
    return readLabel(IStringStream(internalField().group())()) - 1;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sectionalMixedAbsorbingFvPatchScalarField::
sectionalMixedAbsorbingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedAbsorbingFvPatchScalarField(p, iF)
{}

Foam::sectionalMixedAbsorbingFvPatchScalarField::
sectionalMixedAbsorbingFvPatchScalarField
(
    const sectionalMixedAbsorbingFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedAbsorbingFvPatchScalarField(ptf, p, iF, mapper)
{}

Foam::sectionalMixedAbsorbingFvPatchScalarField::
sectionalMixedAbsorbingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedAbsorbingFvPatchScalarField(p, iF, dict)
{}

Foam::sectionalMixedAbsorbingFvPatchScalarField::
sectionalMixedAbsorbingFvPatchScalarField
(
    const sectionalMixedAbsorbingFvPatchScalarField& tppsf
)
:
    mixedAbsorbingFvPatchScalarField(tppsf)
{}

Foam::sectionalMixedAbsorbingFvPatchScalarField::
sectionalMixedAbsorbingFvPatchScalarField
(
    const sectionalMixedAbsorbingFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedAbsorbingFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sectionalMixedAbsorbingFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Set the drift velocity and diffusivity names corresponding to those
    // belonging to this section

    if (db().foundObject<fixedSectionalSystem>("fixedSectionalSystem"))
    {
        const fixedSectionalSystem& system =
            db().lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

        const label sec = this->sectionIndex();
        const word sectionName(system.distribution()[sec].sectionName());

        const word driftVelocityName =
            system.aerosol().drift().dispInertialDrift()
           .velocityFieldName(sectionName);

        // It could be that the drift velocity field is not yet generated at
        // this point. Check this.

        if (db().foundObject<volVectorField>(driftVelocityName))
        {
            this->driftVelocityName_ = driftVelocityName;
        }

        this->diffusivityName_ = system.distribution()[sec].D().name();
    }

    mixedAbsorbingFvPatchScalarField::updateCoeffs();
}

void Foam::sectionalMixedAbsorbingFvPatchScalarField::write(Ostream& os) const
{
    mixedAbsorbingFvPatchScalarField::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        sectionalMixedAbsorbingFvPatchScalarField
    );
}

// ************************************************************************* //
