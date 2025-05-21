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

#include "WilkeWassiljewaComponentMixing.H"
#include "addToRunTimeSelectionTable.H"
#include "rhoAerosolPhaseThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WilkeWassiljewaComponentMixing, 0);
addToRunTimeSelectionTable(componentMixingModel, WilkeWassiljewaComponentMixing, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WilkeWassiljewaComponentMixing::WilkeWassiljewaComponentMixing
(
    const rhoAerosolPhaseThermo& thermo,
    const fieldType& field
)
:
    componentMixingModel(type(), thermo, field)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

WilkeWassiljewaComponentMixing::~WilkeWassiljewaComponentMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void WilkeWassiljewaComponentMixing::mix(volScalarField& phi)
{
    const fvMesh& mesh = thermo_.mesh();

    const label N(thermo_.species().size());

    PtrList<volScalarField> mu(N);
    PtrList<volScalarField> alpha(N);

    const PtrList<volScalarField>& Y = thermo_.composition().Y();

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    const basicSpecieMixture& comp = thermo_.composition();

    // Set molecular weights

    scalarList W(N, 0.0);

    forAll(W, j)
    {
        W[j] = thermo_.composition().W(j);
    }

    forAll(Y, j)
    {
        // Set individual viscosity fields

        mu.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("mu", thermo_.species()[j]),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimViscosity*dimDensity, 0.0)
            )
        );

        volScalarField& muj = mu[j];

        forAll(muj.internalField(), celli)
        {
            muj[celli] = comp.mu(j, p[celli], T[celli]);
        }

        forAll(muj.boundaryField(), patchi)
        {
            scalarField& mujp = muj.boundaryFieldRef()[patchi];

            const scalarField& Tp = T.boundaryField()[patchi];
            const scalarField& pp = p.boundaryField()[patchi];

            forAll(mujp, facei)
            {
                mujp[facei] = comp.mu(j, pp[facei], Tp[facei]);
            }
        }

        // Set individual conductivity fields if this is conductivity

        if (field_ == conductivity)
        {
            alpha.set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("alpha", thermo_.species()[j]),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh,
                    dimensionedScalar("zero", dimViscosity*dimDensity, 0.0)
                )
            );

            volScalarField& alphaj = alpha[j];

            forAll(alphaj.internalField(), celli)
            {
                alphaj[celli] = comp.alphah(j, p[celli], T[celli]);
            }

            forAll(alphaj.boundaryField(), patchi)
            {
                scalarField& alphajp = alphaj.boundaryFieldRef()[patchi];

                const scalarField& Tp = T.boundaryField()[patchi];
                const scalarField& pp = p.boundaryField()[patchi];

                forAll(alphajp, facei)
                {
                    alphajp[facei] = comp.alphah(j, pp[facei], Tp[facei]);
                }
            }
        }
    }

    // Compute mole fractions

    PtrList<volScalarField> x(N);

    volScalarField sumYOverW(Y[0]/W[0]);

    for (label j = 1; j < N; j++)
    {
        sumYOverW += Y[j]/W[j];
    }

    forAll(x, j)
    {
        x.set(j, new volScalarField(Y[j]/W[j]/sumYOverW));
    }

    // Compute phiij matrix (only set upper matrix triangle)

    PtrList<volScalarField> phiij(N*N);

    forAll(x, j)
    {
        forAll(x, k)
        {
            const label i(j*N+k);

            if (j < k)
            {
                phiij.set
                (
                    i,
                    new volScalarField
                    (
                        sqr(1.0+sqrt(mu[j]/mu[k])*pow(W[j]/W[k],0.25))
                      / sqrt(8.0*(1.0+W[j]/W[k]))
                    )
                );
            }
        }
    }

    PtrList<volScalarField>& psi = field_ == viscosity ? mu : alpha;

    phi *= 0.0;

    forAll(x, j)
    {
        // Diagonal component of denominator

        volScalarField sumxkPhijk(x[j]);

        // Off-diagonal components of denominator

        forAll(x, k)
        {
            if (j != k)
            {
                // Swap j and k if j > k

                const label i(j < k ? j*N+k : k*N+j);

                if (j < k)
                {
                    sumxkPhijk += x[k]*phiij[i];
                }
                else
                {
                    // Turn upper triangle value into lower triangle value

                    sumxkPhijk += x[k]*phiij[i]*mu[k]/mu[j]*W[j]/W[k];
                }
            }
        }

        phi += x[j]*psi[j]/sumxkPhijk;
    }

    phi.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
