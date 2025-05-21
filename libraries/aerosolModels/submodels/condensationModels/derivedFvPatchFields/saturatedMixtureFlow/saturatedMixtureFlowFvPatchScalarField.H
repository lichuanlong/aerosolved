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

\file saturatedMixtureFlowFvPatchScalarField.H

\brief Saturated mixture flow BC

This boundary condition computes the equilibrium mixture that would arrise from
a well mixed continuous and dipersed flow. Given the composition of continuous
phase and dispersed phase, and a quality of flow (defined as the dispersed mass
flow rate over the total mass flow rate), the equilibrium mixture mass fractions
are calculated.

This boundary condition requires all mass fraction fields to have this boundary
condition specified. If the quality of the flow is larger than zero, at least
one specie must be part of the dispersed phase. Also, at least one specie must
be part of continuous phase.

*/

#ifndef saturatedMixtureFlowFvPatchScalarField_H
#define saturatedMixtureFlowFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "aerosolModel.H"
#include "autoPtr.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class saturatedMixtureFlowFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class saturatedMixtureFlowFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

    // Private Data

        //- Iterative solution tolerance, defaults to 1e-8
        scalar tolerance_;

        //- Maximum number of iterations, defaults to 1000
        label maxIter_;

        //- Quality of the flow (dispersed mass flow over total mass flow)
        scalar quality_;

        //- List of continuous species volume fractions
        HashTable<scalar> alphaCont_;

        //- List of dispersed species volume fractions (active species only)
        HashTable<scalar> alphaDisp_;

        //- List of activity coefficients
        HashTable<scalar> gamma_;

        //- Switch if fields are checked
        Switch checkFields_;

    // Private Member Functions

        //- Function to set equilibrium mole fractions at given face. Retun true
        //  if converged.
        bool solve
        (
            PtrList<scalarField>& X,
            PtrList<scalarField>& W,
            const PtrList<scalarField>& S,
            const label facei
        ) const;

        //- Check if all fields have this BC
        void checkFields();

public:

    //- Runtime type information
    TypeName("saturatedMixtureFlow");


    // Constructors

        //- Construct from patch and internal field
        saturatedMixtureFlowFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        saturatedMixtureFlowFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given saturatedMixtureFlowFvPatchScalarField
        // onto a new patch
        saturatedMixtureFlowFvPatchScalarField
        (
            const saturatedMixtureFlowFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        saturatedMixtureFlowFvPatchScalarField
        (
            const saturatedMixtureFlowFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new saturatedMixtureFlowFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        saturatedMixtureFlowFvPatchScalarField
        (
            const saturatedMixtureFlowFvPatchScalarField&,
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
                new saturatedMixtureFlowFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

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
