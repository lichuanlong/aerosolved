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

#include "accumulateFlux.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(accumulateFlux, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        accumulateFlux,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::accumulateFlux::calc()
{
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(fieldName_);

    volScalarField& accumulated = accumulatedPtr_();

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];

        if (!mesh_.boundaryMesh()[patchi].coupled())
        {
            accumulated.boundaryFieldRef()[patchi] +=
                mesh_.time().deltaTValue()
              * phi.boundaryField()[patchi]
              / patch.magFaceAreas();
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::accumulateFlux::accumulateFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "phi")
{
    setResultName(typeName, "phi");

    // Create the accumulation field

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(fieldName_);

    accumulatedPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                resultName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                resultName_,
                phi.dimensions()*dimTime/dimArea,
                0.0
            )
        )
    );
}

bool Foam::functionObjects::accumulateFlux::execute()
{
    return calc();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::accumulateFlux::~accumulateFlux()
{}


// ************************************************************************* //
