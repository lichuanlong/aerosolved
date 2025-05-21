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

#include "PecletNumber.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(PecletNumber, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        PecletNumber,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::PecletNumber::calc()
{
    const volVectorField& V =
        mesh_.lookupObject<volVectorField>(driftVelocityName_);

    const volScalarField& D =
        mesh_.lookupObject<volScalarField>(diffusivityName_);

    const volScalarField W(mag(aerosol_.U() + V));

    const tmp<volScalarField> tPe
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
            dimensionedScalar("Pe", dimless, 0.0)
        )
    );

    volScalarField& Pe = tPe.ref();

    // For cells, approximate the cell dimension from its cube root volume

    Pe.field() =
        pow(mesh_.V().field(), 1.0/3.0)
      * W.field()
      / max(D.field(),SMALL);

    // For faces, use the cell-to-face distance as length scale

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const fvPatch& p = W.boundaryField()[patchi].patch();

        const scalarField& Wp = W.boundaryField()[patchi];
        const scalarField& Dp = D.boundaryField()[patchi];

        Pe.boundaryFieldRef()[patchi] = mag(p.delta())*Wp/max(Dp,SMALL);
    }

    return store(resultName_, tPe, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::PecletNumber::PecletNumber
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    driftVelocityName_(dict.lookupOrDefault<word>("driftVelocity", "Ur")),
    diffusivityName_(dict.lookupOrDefault<word>("diffusivity", "DDisp")),
    aerosol_(lookupObject<aerosolModel>("aerosolProperties")),
    resultName_
    (
        dict.lookupOrDefault<word>
        (
            "result",
            word("Pe(" + driftVelocityName_ + "," + diffusivityName_ + ")")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::PecletNumber::~PecletNumber()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::PecletNumber::execute()
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


bool Foam::functionObjects::PecletNumber::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::PecletNumber::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
