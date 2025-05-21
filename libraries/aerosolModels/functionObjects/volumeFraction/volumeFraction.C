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

#include "volumeFraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeFraction, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        volumeFraction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::volumeFraction::calc()
{
    aerosolThermo& thermo = aerosol_.thermo();

    const volScalarField rho(thermo.rho());

    const label j(thermo.species()[speciesName_]);

    const volScalarField& massFrac =
        phaseName_ == thermo.phaseNameCont()
      ? thermo.Y()[j]
      : thermo.Z()[j];

    volScalarField rhoSpecies(rho*0.0);

    if (phaseName_ == thermo.phaseNameCont())
    {
        rhoSpecies.primitiveFieldRef() =
            thermo.rhoCont(speciesName_);

        forAll(rhoSpecies.boundaryField(), patchi)
        {
            rhoSpecies.boundaryFieldRef()[patchi] =
                thermo.rhoCont(speciesName_, patchi);
        }
    }
    else
    {
        rhoSpecies.primitiveFieldRef() =
            thermo.rhoDisp(speciesName_);

        forAll(rhoSpecies.boundaryField(), patchi)
        {
            rhoSpecies.boundaryFieldRef()[patchi] =
                thermo.rhoDisp(speciesName_, patchi);
        }
    }

    tmp<volScalarField> talpha
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
            massFrac*rho/rhoSpecies
        )
    );

    return store(resultName_, talpha, true);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeFraction::volumeFraction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    aerosol_(lookupObjectRef<aerosolModel>("aerosolProperties")),
    speciesName_(dict.lookup("species")),
    phaseName_(dict.lookup("phase")),
    resultName_
    (
        dict.lookupOrDefault<word>
        (
            "result",
            IOobject::groupName
            (
                "alpha",
                IOobject::groupName(speciesName_, phaseName_)
            )
        )
    )
{
    aerosolThermo& thermo = aerosol_.thermo();

    if
    (
        phaseName_ != thermo.phaseNameCont()
     && phaseName_ != thermo.phaseNameDisp()
    )
    {
        FatalErrorInFunction
            << "Invalid phase specified. Phase should be "
            << thermo.phaseNameCont() << " or "
            << thermo.phaseNameDisp() << "." << endl
            << abort(FatalError);
    }

    calc();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeFraction::~volumeFraction()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volumeFraction::execute()
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


bool Foam::functionObjects::volumeFraction::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::volumeFraction::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
