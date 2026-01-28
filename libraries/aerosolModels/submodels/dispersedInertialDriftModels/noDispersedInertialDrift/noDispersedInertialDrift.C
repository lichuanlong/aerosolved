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

#include "noDispersedInertialDrift.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noDispersedInertialDrift, 0);
addToRunTimeSelectionTable(dispersedInertialDriftModel, noDispersedInertialDrift, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noDispersedInertialDrift::noDispersedInertialDrift
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    dispersedInertialDriftModel(type(), aerosol, dict),
    dummy_
    (
        IOobject
        (
            "V.dummy",
            aerosol.mesh().time().timeName(),
            aerosol.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aerosol.mesh(),
        dimensionedVector("V", dimVelocity, vector::zero)
    )
{
    readBaseFieldIfPresent();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noDispersedInertialDrift::~noDispersedInertialDrift()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volVectorField& noDispersedInertialDrift::V
(
    const volScalarField& d,
    const word sizeName
)
{
    return dummy_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
