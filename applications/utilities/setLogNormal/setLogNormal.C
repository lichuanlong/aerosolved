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

\file setLogNormal.C
\brief Utility to set an M field given a count median diameter (CMD) and
logarithm of the geometric standard deviation (sigma) assuming a log-normal size
distribution

*/


#include "argList.H"
#include "fvCFD.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("sigma");
    argList::validArgs.append("CMD");

    #include "setRootCase.H"

    const dimensionedScalar sigma
    (
        "sigma",
        dimless,
        args.get<scalar>(1)
    );

    const dimensionedScalar CMD
    (
        "CMD",
        dimLength,
        args.get<scalar>(2)
    );

    #include "createTime.H"
    #include "createMesh.H"

    aerosolThermo thermo(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    volScalarField M
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const speciesTable& dispSpecies = thermo.thermoDisp().species();

    const volScalarField rhol(thermo.thermoDisp().rho());

    const scalar pi = constant::mathematical::pi;

    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("alpha", dimless, 0.0)
    );

    forAll(dispSpecies, j)
    {
        alpha += thermo.Z()[j];
    }

    M = 6.0*alpha*Foam::exp(-4.5*Foam::sqr(Foam::log(sigma)))
          / (pi*rhol*Foam::pow(CMD,3.0));

    M.write();

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
