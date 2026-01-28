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

\file mixedAbsorbingFvPatchScalarField.H

\brief Mixed BC for droplet absorption driven by both convection and diffusion

Droplets can deposition on a wall by both convection (drift) and diffusion.
Large particle deposition is driven by drift, for which a homogeneous Neumann BC
is needed. Small particle deposition is driven by diffusion, for which a
homogeneous Dirichlet BC is needed. The mixedAbsorbing BC switches between these
two limits, based on the ratio  \f$r = |V_c|/|V_d|\f$ with \f$V_c\f$ a measure
for the convective velocity and \f$V_d\f$ a measure for the diffusive velocity.
The ratio is then used in the function \f$f = r^\gamma/(r^\gamma+1)\f$, which is
in [0,1] and switches between Dirichlet and Neumann.

By default, Homogeneous Dirichlet is strictly applied when the convective flux
is inward, in order to avoid particles being entrained from the boundary. This
can be disabled by setting the clipInwardFlux switch to false.

*/

#ifndef mixedAbsorbingFvPatchScalarField_H
#define mixedAbsorbingFvPatchScalarField_H

#include "aerosolModel.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class MixedAbsorbinFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mixedAbsorbingFvPatchScalarField
:
    public mixedFvPatchScalarField
{
protected:

    // Protected Data

        //- Name of the dispersed drift velocity field (defaults to Ur)
        word driftVelocityName_;

        //- Name of the diffusivity field (defaults to DDisp)
        word diffusivityName_;

        //- Transition parameter (defaults to unity)
        scalar gamma_;

        //- Clip inward flux (defaults to true)
        bool clipInwardFlux_;


    // Protected Member Functions

        //- Return the Peclet number field
        virtual tmp<scalarField> Peclet() const;

        //- Inward convective flux mask, returning unity if the convective flux
        //  is inward
        virtual tmp<scalarField> inwardFluxMask() const;


public:

    //- Runtime type information
    TypeName("mixedAbsorbing");

    // Constructors

        //- Construct from patch and internal field
        mixedAbsorbingFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixedAbsorbingFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedAbsorbingFvPatchScalarField
        // onto a new patch
        mixedAbsorbingFvPatchScalarField
        (
            const mixedAbsorbingFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixedAbsorbingFvPatchScalarField
        (
            const mixedAbsorbingFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mixedAbsorbingFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mixedAbsorbingFvPatchScalarField
        (
            const mixedAbsorbingFvPatchScalarField&,
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
                new mixedAbsorbingFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

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
