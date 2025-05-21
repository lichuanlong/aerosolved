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

\file subGridDepositionVelocityFvPatchVectorField.H \brief Abstract sub-grid
deposition velocity BC

Uses the subGridDepositionModel to determine the deposition velocity at the
patch. Classes that inherit this BC should specify the diameter using the d()
function

*/

#ifndef subGridDepositionVelocityFvPatchVectorField_H
#define subGridDepositionVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class subGridDepositionVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class subGridDepositionVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Gravitational acceleration vector
        vector g_;

        //- Maximum number of Newton iterations
        label maxIter_;

        //- Newton iteration tolerance
        scalar tolerance_;


public:

    //- Runtime type information
    TypeName("subGridDepositionVelocity");

    // Constructors

        //- Construct from patch and internal field
        subGridDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        subGridDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch  onto a new patch
        subGridDepositionVelocityFvPatchVectorField
        (
            const subGridDepositionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        subGridDepositionVelocityFvPatchVectorField
        (
            const subGridDepositionVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );


    // Member functions

        //- Calculate the mean particle diameter at the patch
        virtual tmp<scalarField> d() const = 0;

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
