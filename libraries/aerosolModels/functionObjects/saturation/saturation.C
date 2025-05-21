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

#include "saturation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(saturation, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        saturation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::saturation::calc()
{
    aerosolThermo& thermo = aerosol_.thermo();

    const scalarField sumY(thermo.sumY().ref().field());

    PtrList<scalarField> y(thermo.species().size());

    forAll(y, j)
    {
        y.set(j, new scalarField(thermo.Y()[j].field()/sumY));
    }

    scalarList W(thermo.species().size(), 0.0);

    forAll(W, j)
    {
        W[j] = thermo.thermoCont().composition().W(j);
    }

    scalarField sumyOverW(sumY.size(), 0.0);

    forAll(y, j)
    {
        sumyOverW += y[j]/W[j];
    }

    const scalarField pSat(thermo.pSat(speciesName_));

    const scalarField& p = thermo.p().field();

    tmp<volScalarField> tS
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
            dimensionedScalar("S", dimless, 0.0)
        )
    );

    const label j(thermo.species()[speciesName_]);

    tS.ref().field() = p*y[j]/W[j]/sumyOverW/pSat;

    return store(resultName_, tS, true);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::saturation::saturation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    aerosol_(lookupObjectRef<aerosolModel>("aerosolProperties")),
    speciesName_(dict.lookup("species")),
    resultName_
    (
        dict.lookupOrDefault<word>
        (
            "result",
            IOobject::groupName("S", speciesName_)
        )
    )
{
    calc();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::saturation::~saturation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::saturation::execute()
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


bool Foam::functionObjects::saturation::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::saturation::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
