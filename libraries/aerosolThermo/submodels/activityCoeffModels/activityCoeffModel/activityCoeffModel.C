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

#include "activityCoeffModel.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(activityCoeffModel, 0);
defineRunTimeSelectionTable(activityCoeffModel, dictionary);

const word activityCoeffModel::dictName = "activityCoeff";

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

PtrList<scalarField> activityCoeffModel::unityGamma
(
    const label nSpecies,
    const label nCells
) const
{
    PtrList<scalarField> gamma(nSpecies);

    forAll(gamma, i)
    {
        gamma.set(i, new scalarField(nCells, 1.0));
    }

    return gamma;
}

tmp<scalarField> activityCoeffModel::sumList
(
    const PtrList<scalarField>& x
) const
{
    tmp<scalarField> tSum(new scalarField(x[0]));

    for (int i = 1; i < x.size(); i++)
    {
        tSum.ref() += x[i];
    }

    return tSum;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

activityCoeffModel::activityCoeffModel(aerosolThermo& thermo)
:
    thermo_(thermo),
    dict_
    (
        thermo.found(activityCoeffModel::dictName)
      ? thermo.subDict(activityCoeffModel::dictName)
      : dictionary::null
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

activityCoeffModel::~activityCoeffModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

PtrList<scalarField> activityCoeffModel::gamma
(
    const speciesTable& species
) const
{
    // Compute the molar liquid composition

    const basicSpecieMixture& compCont = thermo_.thermoCont().composition();
    const speciesTable& activeSpecies = thermo_.activeSpecies();

    const scalarField sumZ(min(thermo_.sumZ()->field(),1.0));

    scalarList W(activeSpecies.size());

    forAll(activeSpecies, j)
    {
        W[j] = compCont.W(j);
    }

    PtrList<scalarField> x(activeSpecies.size());

    forAll(x, i)
    {
        x.set(i, new scalarField(sumZ.size(), 0.0));
    }

    forAll(sumZ, celli)
    {
        if (sumZ[celli] > VSMALL)
        {
            scalarList Z(activeSpecies.size());

            forAll(activeSpecies, j)
            {
                Z[j] = thermo_.Z()[j][celli];
            }

            const scalarList z(Z/sumZ[celli]);

            forAll(activeSpecies, j)
            {
                x[j][celli] = z[j]/W[j]/sum(z/W);
            }
        }
    }

    return
        gamma
        (
            species,
            thermo_.p().field(),
            thermo_.T().field(),
            x
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
