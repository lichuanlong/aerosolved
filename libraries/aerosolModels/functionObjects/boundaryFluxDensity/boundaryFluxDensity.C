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

#include "boundaryFluxDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(boundaryFluxDensity, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        boundaryFluxDensity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::boundaryFluxDensity::calc()
{
    const surfaceScalarField& phi =
        lookupObject<surfaceScalarField>(phiName_);

    const tmp<volScalarField> tflux
    (
        new volScalarField
        (
            IOobject
            (
                "t" + resultName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", phi.dimensions()/dimArea, 0)
        )
    );

    volScalarField& flux = tflux.ref();

    forAll(mesh_.boundaryMesh(), patchi)
    {
        flux.boundaryFieldRef()[patchi] =
            phi.boundaryField()[patchi]
          / mesh_.magSf().boundaryField()[patchi];
    }

    return store(resultName_, tflux, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::boundaryFluxDensity::boundaryFluxDensity
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phiName_(dict.lookup("phi")),
    resultName_
    (
        dict.lookupOrDefault<word>
        (
            "result",
            word("boundaryFluxDensity("+phiName_+")")
        )
    )
{
    calc();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::boundaryFluxDensity::~boundaryFluxDensity()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::boundaryFluxDensity::execute()
{
    if (!calc())
    {
        Warning
            << "    functionObjects::" << type() << " " << name()
            << " failed to execute." << endl;

        // Clear the result field from the objectRegistry if present
        clear();

        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::functionObjects::boundaryFluxDensity::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::boundaryFluxDensity::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
