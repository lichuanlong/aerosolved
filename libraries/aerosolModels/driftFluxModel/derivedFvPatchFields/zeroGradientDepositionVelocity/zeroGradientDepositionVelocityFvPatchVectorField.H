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

\file zeroGradientDepositionVelocityFvPatchVectorField.H

\brief BC for the relative velocity, assuring deposition of dispersed phase

This BC applies to the relative velocity V. It sets either the relative velocity
or absolute dispersed phase velocity Ud cell-centered value on the patch. If
this gives a dispersed phase velocity into the domain, it is clipped in such a
way that the dispersed phase velocity becomes zero. If it gives a dispersed
phase velocity that is into the wall, then it is left as such, resulting
effectively in a zero-gradient BC for either V or Ud. The relation between
dispersed phase velocity Ud, mixture velocity U and relative velocity V is given
by Ud = U + V - VCorr, with VCorr the corrective velocity giving the difference
between mixture velocity and continuous phase velocity. By default, the BC
applies the zero-gradient mechanism for the relative velocity, but by setting
the dispersedVelocity switch to true, it applies it to the dispersed phase
velocity.

*/

#ifndef zeroGradientDepositionVelocityFvPatchVectorField_H
#define zeroGradientDepositionVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
        Class zeroGradientDepositionVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class zeroGradientDepositionVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:

    // Private Data

        //- Zero-gradient mode: apply zero-gradient to the relative velocity (V)
        //  or to the dispersed velocity (Ud). Defaults to V.
        word zeroGradientMode_;

        //- Clip inward flux (defaults to true)
        bool clipInwardFlux_;

public:

   //- Runtime type information
   TypeName("zeroGradientDepositionVelocity");


   // Constructors

        //- Construct from patch and internal field
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const zeroGradientDepositionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const zeroGradientDepositionVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new zeroGradientDepositionVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const zeroGradientDepositionVelocityFvPatchVectorField&,
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
                new zeroGradientDepositionVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

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
