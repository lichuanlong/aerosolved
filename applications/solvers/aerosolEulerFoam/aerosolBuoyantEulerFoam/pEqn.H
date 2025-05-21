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

rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
      + rhorAUf*fvc::ddtCorr(rho, U, phi)
    )
  + phig
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

fvScalarMatrix p_rghDDtEqn
(
    fvm::ddt(psi, p_rgh)
  + fvc::ddt(psi, rho)*gh
  + fvc::div(phiHbyA)
  ==
    fvOptions(psi, p_rgh, rho.name())
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        p_rghDDtEqn
      - fvm::laplacian(rhorAUf, p_rgh)
    );

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + p_rghEqn.flux();

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        K = 0.5*magSqr(U);

        aerosol->drift().correctPhaseVelocities();
    }
}

p = p_rgh + rho*gh;

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"
