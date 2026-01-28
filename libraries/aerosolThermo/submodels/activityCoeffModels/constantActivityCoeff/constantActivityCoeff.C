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

#include "constantActivityCoeff.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolThermo.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantActivityCoeff, 0);
addToRunTimeSelectionTable
(
    activityCoeffModel,
    constantActivityCoeff,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantActivityCoeff::constantActivityCoeff(aerosolThermo& thermo)
:
    activityCoeffModel(thermo),
    coeffs_(thermo.activeSpecies().size(), 1.0)
{
    forAll(thermo.activeSpecies(), j)
    {
        coeffs_[j] = dict_.lookupOrDefault<scalar>
            (
                thermo.activeSpecies()[j],
                1.0
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantActivityCoeff::~constantActivityCoeff()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

PtrList<scalarField> constantActivityCoeff::gamma
(
    const speciesTable& species
) const
{
    PtrList<scalarField> gamma
    (
        unityGamma(species.size(), thermo_.mesh().nCells())
    );

    forAll(gamma, i)
    {
        gamma[i] *= coeffs_[thermo_.activeSpecies()[species[i]]];
    }

    return gamma;
}

PtrList<scalarField> constantActivityCoeff::gamma
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
