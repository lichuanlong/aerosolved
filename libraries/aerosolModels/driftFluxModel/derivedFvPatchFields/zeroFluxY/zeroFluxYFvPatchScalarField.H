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

\file zeroFluxYFvPatchScalarField.H
\brief The ZG absorbing wall BC

This BC uses the zero-gradient approach when evaluating the value at the patch,
or the zero-value when evaluating the gradient at the patch. This is useful to
allow for both convection and diffusion of the scalar at hand.

*/

#ifndef zeroFluxYFvPatchScalarFields_H
#define zeroFluxYFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class zeroFluxYFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class zeroFluxYFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("zeroFluxY");


    // Constructors


        //- Construct from patch and internal field
        zeroFluxYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroFluxYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given zeroFluxYFvPatchScalarField
        // onto a new patch
        zeroFluxYFvPatchScalarField
        (
            const zeroFluxYFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        zeroFluxYFvPatchScalarField
        (
            const zeroFluxYFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new zeroFluxYFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroFluxYFvPatchScalarField
        (
            const zeroFluxYFvPatchScalarField&,
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
                new zeroFluxYFvPatchScalarField(*this, iF)
            );
        }



    // Member functions

        // Evaluation functions

        // Get the section index for the Z
        Foam::label speciesIndex() const;

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
