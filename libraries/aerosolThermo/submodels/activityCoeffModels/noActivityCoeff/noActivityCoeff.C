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

#include "noActivityCoeff.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolThermo.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noActivityCoeff, 0);
addToRunTimeSelectionTable
(
    activityCoeffModel,
    noActivityCoeff,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noActivityCoeff::noActivityCoeff(aerosolThermo& thermo)
:
    activityCoeffModel(thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noActivityCoeff::~noActivityCoeff()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

PtrList<scalarField> noActivityCoeff::gamma
(
    const speciesTable& species
) const
{
    return unityGamma(species.size(), thermo_.mesh().nCells());
}

PtrList<scalarField> noActivityCoeff::gamma
(
    const speciesTable& species,
    const scalarField&,
    const scalarField&,
    const PtrList<scalarField>&
) const
{
    return gamma(species);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
