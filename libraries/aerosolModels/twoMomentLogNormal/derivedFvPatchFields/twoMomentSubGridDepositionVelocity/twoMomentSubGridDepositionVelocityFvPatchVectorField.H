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

\file twoMomentSubGridDepositionVelocityFvPatchVectorField.H
\brief Sectional sub-grid deposition velocity BC

Uses the subGridDepositionModel to determine the deposition velocity at the
patch.

*/

#ifndef twoMomentSubGridDepositionVelocityFvPatchVectorField_H
#define twoMomentSubGridDepositionVelocityFvPatchVectorField_H

#include "subGridDepositionVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class twoMomentSubGridDepositionVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class twoMomentSubGridDepositionVelocityFvPatchVectorField
:
    public subGridDepositionVelocityFvPatchVectorField
{
    // Private member functions

        //- Find the section index to which this BC's field belongs to
        label sectionIndex() const;


public:

    //- Runtime type information
    TypeName("twoMomentSubGridDepositionVelocity");


    // Constructors

        //- Construct from patch and internal field
        twoMomentSubGridDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        twoMomentSubGridDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch  onto a new patch
        twoMomentSubGridDepositionVelocityFvPatchVectorField
        (
            const twoMomentSubGridDepositionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new twoMomentSubGridDepositionVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        twoMomentSubGridDepositionVelocityFvPatchVectorField
        (
            const twoMomentSubGridDepositionVelocityFvPatchVectorField&,
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
                new twoMomentSubGridDepositionVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Calculate the mean particle diameter at the patch from the section
        virtual tmp<scalarField> d() const;

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
