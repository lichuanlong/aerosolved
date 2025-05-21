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

/**

\file sectionalMixedAbsorbingFvPatchScalarField.H

\brief Mixed BC for sectional absorption driven by both convection and diffusion

This BC is based on the mixedAbsorbing BC. All the sectionalMixedAbsorbing BC
does is overwrite phiName and diffusivityName with the appropriate names
belonging to the section to which this BC is applied.

*/

#ifndef sectionalMixedAbsorbingFvPatchScalarField_H
#define sectionalMixedAbsorbingFvPatchScalarField_H

#include "aerosolModel.H"
#include "mixedAbsorbingFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class MixedAbsorbinFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sectionalMixedAbsorbingFvPatchScalarField
:
    public mixedAbsorbingFvPatchScalarField
{
private:

    // Private Member Functions

        //- Index of this section
        label sectionIndex() const;

public:

    //- Runtime type information
    TypeName("sectionalMixedAbsorbing");

    // Constructors

        //- Construct from patch and internal field
        sectionalMixedAbsorbingFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sectionalMixedAbsorbingFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given sectionalMixedAbsorbingFvPatchScalarField
        // onto a new patch
        sectionalMixedAbsorbingFvPatchScalarField
        (
            const sectionalMixedAbsorbingFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sectionalMixedAbsorbingFvPatchScalarField
        (
            const sectionalMixedAbsorbingFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new sectionalMixedAbsorbingFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sectionalMixedAbsorbingFvPatchScalarField
        (
            const sectionalMixedAbsorbingFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new sectionalMixedAbsorbingFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
