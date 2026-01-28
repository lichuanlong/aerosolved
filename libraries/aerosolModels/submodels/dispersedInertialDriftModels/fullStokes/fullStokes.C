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

#include "fullStokes.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fullStokes, 0);
addToRunTimeSelectionTable(dispersedInertialDriftModel, fullStokes, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fullStokes::fullStokes
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    dispersedInertialDriftModel(type(), aerosol, dict),
    Y_
    (
        IOobject
        (
            "Y",
            aerosol.mesh().time().timeName(),
            aerosol.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aerosol.thermo().sumY()
    ),
    maxIter_(readScalar(dict.lookup("maxIter"))),
    TOL_(readScalar(dict.lookup("tolerance")))
{
    this->readBaseField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fullStokes::~fullStokes()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volVectorField& fullStokes::V
(
    const volScalarField& d,
    const word sizeName
)
{
    volVectorField& V = this->getCachedField(sizeName);

    const fvMesh& mesh = aerosol_.mesh();

    const aerosolThermo& thermo = aerosol_.thermo();

    const volScalarField& rho = aerosol_.rho();
    const volScalarField& rhol = thermo.thermoDisp().rho();
    const volScalarField& rhog = thermo.thermoCont().rho();

    const volVectorField& U = aerosol_.U();
    const surfaceScalarField& phi = aerosol_.phi();

    const volScalarField& mug = thermo.thermoCont().mu();

    const volScalarField tau(Foam::sqr(d)*rhol/(18.0*mug));
    const volVectorField G((rhol-rhog)/rhol*g_);

    const surfaceScalarField phiU(phi/linearInterpolate(rho));

    Y_ = thermo.sumY();

    const int lvl = Info.level;

    scalar r0(0.0);

    for (label iter = 0; iter < maxIter_; iter++)
    {
        const surfaceScalarField phi(phiU + fvc::flux(Y_*V));
        // Schiller And Naumann (1935) Ver. Deutsch Ing.
        const volScalarField f(1.0 + 0.15 * pow(Re(d,V),0.687)  );

        fvVectorMatrix VEqn
        (
            fvm::ddt(Y_,V)
          + fvc::ddt(U)

          + fv::gaussConvectionScheme<vector>
            (
                mesh,
                phi,
                upwind<vector>(mesh, phi)
            ).fvmDiv(phi*fvc::interpolate(Y_), V)

            // makes sure that this uses the same discretization as the
            // momentum equation. Otherwise: tears and oscillations
          + fvc::div(phi,U,"div(phi,U)")

          - fvm::Sp(fvc::div(phi)*Y_, V)
          - fvc::Sp(fvc::div(phi), U)
          ==
          - fvm::Sp(f/tau, V)
          + G
        );

        VEqn.relax();

        Info.level = 0;

        const vector rv(VEqn.solve(mesh.solver("V")).initialResidual());

        const scalar r(max(max(rv[0], rv[1]), rv[2]));

        if (iter == 0)
        {
            r0 = r;
        }

        Info.level = lvl;

        if (r < TOL_ || iter == (maxIter_-1))
        {
            const scalar maxRe(gMax(Re(d, V)().field()));

            Info<< "fullStokes: Solving for " << V.name()
                << ", max(Re) = " << maxRe << ", Initial residual = " << r0
                << ", Final residual = " << r
                << ", No Iterations " << iter+1 << endl;

            break;
        }
    }

    limit(V);

    return V;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
