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

#include "sectionalPecletNumbers.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedSectionalSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sectionalPecletNumbers, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sectionalPecletNumbers,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::sectionalPecletNumbers::calc()
{
    const fixedSectionalSystem& system =
        mesh_.lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

    bool ret = true;

    forAll(system.distribution(), i)
    {
        const word sectionName(system.distribution()[i].sectionName());

        this->driftVelocityName_ =
            system.aerosol().drift().dispInertialDrift()
           .velocityFieldName(sectionName);

        this->diffusivityName_ = system.distribution()[i].D().name();

        this->resultName_ = IOobject::groupName("Pe", sectionName);

        // Use the calc() function of the PecletNumber class to calculate and
        // store the Peclet number

        ret = PecletNumber::calc() && ret;
    }

    return ret;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalPecletNumbers::sectionalPecletNumbers
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    PecletNumber(name, runTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalPecletNumbers::~sectionalPecletNumbers()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sectionalPecletNumbers::write()
{
    const fixedSectionalSystem& system =
        mesh_.lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

    bool ret = true;

    forAll(system.distribution(), i)
    {
        const word sectionName(system.distribution()[i].sectionName());

        this->resultName_ = IOobject::groupName("Pe", sectionName);

        ret = PecletNumber::write() && ret;
    }

    return ret;
}


bool Foam::functionObjects::sectionalPecletNumbers::clear()
{
    const fixedSectionalSystem& system =
        mesh_.lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

    bool ret = true;

    forAll(system.distribution(), i)
    {
        const word sectionName(system.distribution()[i].sectionName());

        this->resultName_ = IOobject::groupName("Pe", sectionName);

        ret = PecletNumber::clear() && ret;
    }

    return ret;
}

// ************************************************************************* //
