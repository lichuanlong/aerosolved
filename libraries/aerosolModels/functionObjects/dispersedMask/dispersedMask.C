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

#include "dispersedMask.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dispersedMask, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        dispersedMask,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::dispersedMask::calc()
{
    const tmp<volScalarField> tdispersedMask
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
            pos(aerosol_.thermo().sumZ()-threshold_)
        )
    );

    return store(resultName_, tdispersedMask, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dispersedMask::dispersedMask
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    threshold_(dict.lookupOrDefault<scalar>("threshold", 1e-6)),
    aerosol_(lookupObject<aerosolModel>("aerosolProperties")),
    resultName_
    (
        dict.lookupOrDefault<word>("result", "dispersedMask")
    )
{
    calc();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::dispersedMask::~dispersedMask()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dispersedMask::execute()
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


bool Foam::functionObjects::dispersedMask::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::dispersedMask::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
