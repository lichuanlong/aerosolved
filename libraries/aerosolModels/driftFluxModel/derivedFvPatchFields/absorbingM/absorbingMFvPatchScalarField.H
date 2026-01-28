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

\file absorbingMFvPatchScalarField.H
\brief Constant sectional BC

*/

#ifndef absorbingMFvPatchScalarField_H
#define absorbingMFvPatchScalarField_H

#include "aerosolModel.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class absorbingMFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class absorbingMFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
 	// Boundary Condition Parametrs

	// Rebound Coefficient
	scalar alpha_= 1.0;

	// Stickiness coefficient
	scalar beta_ = 1.0;

        // Clip M to 0 if inward velocity at patch
        bool clipInwardFlux_ = true;

public:

	//- Runtime type information
	TypeName("absorbingM");


	// Constructors

        //- Construct from patch and internal field
        absorbingMFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        absorbingMFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given absorbingMFvPatchScalarField
        // onto a new patch
        absorbingMFvPatchScalarField
        (
            const absorbingMFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        absorbingMFvPatchScalarField
        (
            const absorbingMFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new absorbingMFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        absorbingMFvPatchScalarField
        (
            const absorbingMFvPatchScalarField&,
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
                new absorbingMFvPatchScalarField(*this, iF)
            );
        }


        // Member functions

        // If fixedSectional, get section index of the droplet
        label  sectionIndex()  const;

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
