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

scalar aerosolCoNum = 0.0;
scalar driftCoNum = 0.0;
scalar relCoNum = 0.0;
scalar corrCoNum = 0.0;

{
    // Aerosol Courant number

    forAll(Z, j)
    {
        aerosolCoNum = max
        (
            gMax
            (
                mag
                (
                    aerosol->R(Z[j])().source()
                  / (rho.field()*mesh.V().field())
                  * runTime.deltaTValue()
                )
            ),
            aerosolCoNum
        );
    }

    const surfaceScalarField& phiInertial = aerosol->drift().phiInertial();
    const surfaceScalarField& phiCorr = aerosol->drift().phiCorr();

    // Drift Courant number

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi+phiInertial))().primitiveField()
      / rho.primitiveField()
    );

    driftCoNum =  0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    // Relative Courant number

    scalarField sumPhiv
    (
        fvc::surfaceSum(mag(phiInertial))().primitiveField()
      / rho.primitiveField()
    );

    relCoNum = 0.5*gMax(sumPhiv/mesh.V().field())*runTime.deltaTValue();

    // Corrective flux Courant number

    scalarField sumPhic
    (
        fvc::surfaceSum(mag(phiCorr))().primitiveField()
      / rho.primitiveField()
    );

    corrCoNum = 0.5*gMax(sumPhic/mesh.V().field())*runTime.deltaTValue();
}

// Print CFL numbers

Info<< "Aerosol Courant Number max: " << aerosolCoNum << nl;
Info<< "Drift Courant Number max: " << driftCoNum << nl;
Info<< "Relative Courant Number max: " << relCoNum << nl;
Info<< "Corrective Courant Number max: " << corrCoNum << nl;

Info<< endl;

// ************************************************************************* //
