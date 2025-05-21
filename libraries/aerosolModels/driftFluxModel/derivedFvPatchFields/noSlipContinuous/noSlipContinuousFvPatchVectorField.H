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

\file noSlipContinuousFvPatchVectorField.H

\brief Mixture velocity boundary condition so that the continuous phase is
no-slip

The mixture velocity is defined as the velocity of average mass. By setting it
equal to dispersed phase velocity, it is implied that the continuous phase
velocity is zero.

*/

#ifndef noSlipContinuousFvPatchVectorField_H
#define noSlipContinuousFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class noSlipContinuousFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class noSlipContinuousFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("noSlipContinuous");


    // Constructors

        //- Construct from patch and internal field
        noSlipContinuousFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        noSlipContinuousFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given noSlipContinuousFvPatchVectorField
        //  onto a new patch
        noSlipContinuousFvPatchVectorField
        (
            const noSlipContinuousFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        noSlipContinuousFvPatchVectorField
        (
            const noSlipContinuousFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new noSlipContinuousFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        noSlipContinuousFvPatchVectorField
        (
            const noSlipContinuousFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new noSlipContinuousFvPatchVectorField(*this, iF)
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
