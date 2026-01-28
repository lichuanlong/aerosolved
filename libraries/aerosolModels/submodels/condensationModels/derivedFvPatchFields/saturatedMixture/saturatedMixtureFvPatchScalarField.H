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

\file saturatedMixtureFvPatchScalarField.H

\brief Saturated mixture BC

This BC can be used to specify a saturated mixture at a patch. We have a set A
of active species having this BC, a set I of inactive species having this BC and
a set O of other species having another BC. The saturation is computed using a
hypothetical pool of liquid, of which the composition determines the continuous
mass fractions of species in A and I. The pool composition for O-species follows
from the evaluation of their respective BCs and Raoult's law. The pool
composition for A-species follows then from the remainder, where the composition
of A-species mixture in the pool liquid is set using the f parameter in the BC.
If there is only one specie in A, then f defaults to unity. From the pool
composition, the vapor mass fractions can be computed using Raoult's law. The
user can specify a total saturation S to have an over- or under-saturated
mixture. S defaults to unity. If there are more than one species in I, the user
should specify the composition of I-species in the continuous mixture of
I-species, using the g parameter. If there is only one specie in I, g defaults
to unity.

*/

#ifndef saturatedMixtureFvPatchScalarField_H
#define saturatedMixtureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "aerosolModel.H"
#include "autoPtr.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class saturatedMixtureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class saturatedMixtureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

    // Private Data

        //- Total saturation fraction
        scalar S_;

        //- Iterative solution tolerance, defaults to 1e-8
        scalar tolerance_;

        //- Switch if fields are checked
        Switch checkFields_;

        //- Alternative saturation temperature (defaults to -1, meaning that the
        //  patch temperature will be taken)
        scalar T_;

        //- List of pool phase active specie mixture mole fractions read
        //  from the f field in the BC dict. It defaults to unity if there is
        //  only one 'type 0' specie
        HashTable<scalar> f_;

        //- List of continuous phase inactive specie mixture mole fractions read
        //  from the g field in the BC dict. It defaults to unity if there is
        //  only one 'type 1' specie
        HashTable<scalar> g_;

        //- List indicating the role each species plays:
        //  0 -> active species having the saturatedMixture BC
        //  1 -> inactive species having the saturatedMixture BC
        //  2 -> other specie not having the saturatedMixture BC
        List<label> specieType_;

        //- Number of species per type
        label N0_;
        label N1_;
        label N2_;


    // Private Member Functions

        //- Check if all fields have this BC
        void checkFields();

public:

    //- Runtime type information
    TypeName("saturatedMixture");


    // Constructors

        //- Construct from patch and internal field
        saturatedMixtureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        saturatedMixtureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given saturatedMixtureFvPatchScalarField
        // onto a new patch
        saturatedMixtureFvPatchScalarField
        (
            const saturatedMixtureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        saturatedMixtureFvPatchScalarField
        (
            const saturatedMixtureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new saturatedMixtureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        saturatedMixtureFvPatchScalarField
        (
            const saturatedMixtureFvPatchScalarField&,
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
                new saturatedMixtureFvPatchScalarField(*this, iF)
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
