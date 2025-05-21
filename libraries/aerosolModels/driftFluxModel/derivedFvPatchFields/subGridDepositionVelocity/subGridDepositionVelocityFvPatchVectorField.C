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

#include "subGridDepositionVelocityFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "aerosolModel.H"
#include "subGridDepositionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subGridDepositionVelocityFvPatchVectorField::
subGridDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    g_(vector::zero),
    maxIter_(99),
    tolerance_(1e-8)
{}


Foam::subGridDepositionVelocityFvPatchVectorField::
subGridDepositionVelocityFvPatchVectorField
(
    const subGridDepositionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    g_(ptf.g_),
    maxIter_(ptf.maxIter_),
    tolerance_(ptf.tolerance_)
{}


Foam::subGridDepositionVelocityFvPatchVectorField::
subGridDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    g_
    (
        uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                patch().boundaryMesh().mesh().time().constant(),
                patch().boundaryMesh().mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedVector("g", dimVelocity/dimTime, vector::zero)
        ).value()
    ),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 99)),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-8))
{
    fvPatchVectorField::operator=(patchInternalField());
}


Foam::subGridDepositionVelocityFvPatchVectorField::
subGridDepositionVelocityFvPatchVectorField
(
    const subGridDepositionVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    g_(fcvpvf.g_),
    maxIter_(fcvpvf.maxIter_),
    tolerance_(fcvpvf.tolerance_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::subGridDepositionVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    fixedValueFvPatchVectorField::evaluate();

    const aerosolModel& aerosol =
        db().lookupObject<aerosolModel>("aerosolProperties");

    const vectorField U
    (
        aerosol.U().boundaryField()[patch().index()].patchInternalField()
    );

    const aerosolThermo& thermo = aerosol.thermo();

    const scalarField rhoc(thermo.thermoCont().rho(patch().index()));
    const scalarField rhod(thermo.thermoDisp().rho(patch().index()));
    const scalarField gamma(rhoc/rhod);

    const scalarField mu(thermo.thermoCont().mu(patch().index()));
    const scalarField d(this->d());
    const scalarField tau(rhod*sqr(d)/(18.0*mu));

    const scalarField delta(mag(patch().delta()));
    const vectorField n(-patch().nf());

    // Solve in terms of absolute velocity W

    vectorField V(patchInternalField());
    vectorField W(V + U);

    forAll(*this, facei)
    {
        // Non-dimensional scalar quantities, projected onto the inward face
        // normal

        const scalar S = delta[facei]/tau[facei];
        const scalar w = (W[facei] & n[facei])/S;

        // Correct the velocity only if the cell-centered velocity is pointing
        // towards the wall and if the particle is not yet touching the wall

        if (w < 0.0 && delta[facei] > d[facei]/2.0)
        {
            // Non-dimensional mixture velocity and gravity

            const scalar u = (U[facei] & n[facei])/S;
            const scalar g =
                (g_ & n[facei]) * (1.0-gamma[facei]) * tau[facei]/S;

            // Construct and solve sub-grid model

            const subGridDepositionModel model
            (
                u,
                g,
                w,
                d[facei]/2.0,
                maxIter_,
                tolerance_
            );

            const scalar wNew =
                model.collision() ? model.v(model.t()) : 0.0;

            // Add correction and translate back to drift velocity

            V[facei] = W[facei] + (wNew - w)*S*n[facei] - U[facei];
        }
    }

    operator==(V);
}


void Foam::subGridDepositionVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);

    os.writeEntryIfDifferent<label>("maxIter", 99, maxIter_);
    os.writeEntryIfDifferent<scalar>("tolerance", 1e-8, tolerance_);

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        subGridDepositionVelocityFvPatchVectorField,
        0
    );
}

// ************************************************************************* //