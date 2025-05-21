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

{
    // Optionally rescale all mass fractions to force a unity sum

    const Switch rescale =
        pimple.dict().lookupOrDefault<Switch>("rescale", true);

    // Optionally add continuous phase diffusion correction to make the
    // diffusion term proportional to mole fraction, not mass fraction

    const Switch molarDiffusion =
        pimple.dict().lookupOrDefault<Switch>("molarDiffusion", false);

    // Prepare fields

    scalarField Yt(mesh.nCells(), 0.0);

    const volScalarField mut(turbulence->mut());
    const volScalarField sumY(thermo.sumY());

    const surfaceScalarField& phiInertial = aerosol->drift().phiInertial();
    surfaceScalarField& phiCorr = aerosol->drift().phiCorr();

    const fv::convectionScheme<scalar>& mvPhi = aerosol->drift().mvPhi();
    const fv::convectionScheme<scalar>& mvPhiCorr =
        aerosol->drift().mvPhiCorr();
    const fv::convectionScheme<scalar>& mvPhiInertial =
        aerosol->drift().mvPhiInertial();

    surfaceScalarField sumPhiEff(0.0*phi);

    // Field to collect all the new drift fluxes in, which is then used to
    // update the phiCorr field

    surfaceScalarField phiCorrNew("phiCorrNew", phiCorr*0.0);

    // Reset the drift stress tensor

    tauDrift *= 0.0;

    // Continuous mass fraction equations

    forAll(Y, j)
    {
        volScalarField& Yj = Y[j];

        const volScalarField& D = aerosol->drift().DCont()[j];

        // Diffusive drift. Add a correction to make diffusion proportional to
        // grad(yj) instead of grad(Yj). This is important for when sum(Zj)
        // becomes significant. Optionally correct the term such that it is
        // proportional to the mole fraction gradient instead of the mass
        // fraction gradient.

        tmp<fvScalarMatrix> tdriftDiff;

        if (!molarDiffusion)
        {
            // Diffusive flux is muEff*grad(yj), yj = phase-specific mass
            // fraction

            surfaceScalarField phiCorrDiff
            (
                linearInterpolate(rho*D)*fvc::snGrad(1.0/sumY)*mesh.magSf()
            );

            tdriftDiff =
            (
              - fvm::laplacian(mut+rho*D/sumY, Yj, "laplacian(mut+rho*D,Yj)")
              - fv::gaussConvectionScheme<scalar>
                (
                    mesh,
                    phiCorrDiff,
                    linear<scalar>(mesh)
                ).fvmDiv(phiCorrDiff, Yj)
            );
        }
        else
        {
            // Diffusive flux is muEff*grad(xj), xj = phase-specific mole
            // fraction

            const volScalarField W(thermo.thermoCont().WMix());

            const dimensionedScalar Wj
            (
                dimMass/dimMoles,
                thermo.thermoCont().composition().W(j)
            );

            const surfaceScalarField phiCorrDiff
            (
                linearInterpolate(rho*D/Wj)*fvc::snGrad(W/sumY)
              * mesh.magSf()
            );

            tdriftDiff =
            (
              - fvm::laplacian
                (
                    mut+rho*D*W/sumY/Wj,
                    Yj,
                    "laplacian(mut+rho*D,Yj)"
                )
              - fv::gaussConvectionScheme<scalar>
                (
                    mesh,
                    phiCorrDiff,
                    linear<scalar>(mesh)
                ).fvmDiv(phiCorrDiff, Yj)
            );
        }

        const fvScalarMatrix& driftDiff = tdriftDiff();

        // Mass fraction equation

        fvScalarMatrix YjEqn
        (
            fvm::ddt(rho, Yj)
          + mvPhi.fvmDiv(phi, Yj)
          + driftDiff
            ==
            aerosol->R(Yj)
          + fvOptions(rho, Yj)
          + mvPhiCorr.fvmDiv(phiCorr, Yj)
        );

        YjEqn.relax();

        fvOptions.constrain(YjEqn);

        YjEqn.solve(mesh.solver("Y"));

        fvOptions.correct(Yj);

        Yj.correctBoundaryConditions();
        Yj.max(0.0);
        Yj.min(1.0);

        Yt += Yj.field();

        // Store fluxes

        phiEffY[j] = YjEqn.flux();
        phiCorrNew += driftDiff.flux();

        sumPhiEff += phiEffY[j];

        // Add to the drift stress tensor, by reconstructing the
        // mixture-relative drift velocity from the fluxes

        tauDrift -=
            sqr
            (
                fvc::reconstruct(driftDiff.flux())
              / max(fvc::average(Yj),1e-9)
              - fvc::reconstruct(phiCorr)
            )
          * Yj/rho;
    }

    // Corrective flux to account for the fact that dispersed phase diffusion is
    // proportional to the gradient of number concentration per unit volume, not
    // number concentration per unit mass

    const volScalarField& D = aerosol->drift().DDisp();
    const surfaceScalarField& phiCorrDiff = aerosol->drift().phiCorrDiff();

    // Dispersed mass fraction equations

    forAll(Z, j)
    {
        volScalarField& Zj = Z[j];

        // Diffusive drift

        fvScalarMatrix driftDiff
        (
          - fvm::laplacian(mut+rho*D, Zj, "laplacian(mut+rho*D,Zj)")
          - fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCorrDiff,
                linear<scalar>(mesh)
            ).fvmDiv(phiCorrDiff, Zj)
        );

        // Inertial drift

        fvScalarMatrix driftInertial
        (
            mvPhiInertial.fvmDiv(phiInertial, Zj)
        );

        // Mass fraction equation

        fvScalarMatrix ZjEqn
        (
            fvm::ddt(rho, Zj)
          + mvPhi.fvmDiv(phi, Zj)
          + driftDiff
          + driftInertial
          ==
            aerosol->R(Zj)
          + fvOptions(rho, Zj)
          + mvPhiCorr.fvmDiv(phiCorr, Zj)
        );

        ZjEqn.relax();

        fvOptions.constrain(ZjEqn);

        ZjEqn.solve(mesh.solver("Y"));

        fvOptions.correct(Zj);

        Zj.max(0.0);
        Zj.min(1.0);
        Zj.correctBoundaryConditions();

        Yt += Zj;

        // Store fluxes

        phiEffZ[j] = ZjEqn.flux();
        phiCorrNew += driftDiff.flux() + driftInertial.flux();

        sumPhiEff += phiEffZ[j];

        // Add to the drift stress tensor, by reconstructing the
        // mixture-relative drift velocity from the fluxes

        tauDrift -=
            sqr
            (
                fvc::reconstruct(driftDiff.flux())
              / max(fvc::average(Zj), 1e-9)
              + fvc::reconstruct(phiInertial-phiCorr)
            )
          * Zj/rho;
    }

    Info<< "Cumulative mass fraction "
        << "min = " << min(Yt) << ", max = " << max(Yt)
        << ", mean defect = " << average(mag(1.0-Yt)) << endl;

    if (rescale)
    {
        Info << "Rescaling mass fractions" << endl;

        forAll(Y, j)
        {
            Y[j].field() /= Yt;
            Y[j].correctBoundaryConditions();
        }

        forAll(Z, j)
        {
            Z[j].field() /= Yt;
            Z[j].correctBoundaryConditions();
        }
    }

    phiCorr = phiCorrNew;
}
