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

#include "Zhang.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolThermo.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Zhang, 0);
addToRunTimeSelectionTable(activityCoeffModel, Zhang, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Zhang::Zhang(aerosolThermo& thermo)
:
    activityCoeffModel(thermo),
    firstSpecieName_(dict_.lookup("firtSpecieName")),
    C1_(readScalar(dict_.lookup("C1"))),
    C2_(readScalar(dict_.lookup("C2")))
{
    if (thermo.activeSpecies().size() != 2)
    {
        FatalErrorInFunction
            << "This activity coefficient model only works for two "
            << "active species" << nl << exit(FatalError);
    }

    if (!thermo.activeSpecies().found(firstSpecieName_))
    {
        FatalErrorInFunction
            << "Could not find specie " << firstSpecieName_
            << " in the list of active species" << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Zhang::~Zhang()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

PtrList<scalarField> Zhang::gamma
(
    const speciesTable& species,
    const scalarField& p,
    const scalarField& T,
    const PtrList<scalarField>& x
) const
{
    if (species.size() != 2)
    {
        FatalErrorInFunction
            << "The " << this->typeName << " activity coefficient model only "
            << "works for two active species."
            << endl << abort(FatalError);
    }

    rhoAerosolPhaseThermo& thermoCont = thermo_.thermoCont();

    const speciesTable& activeSpecies = thermo_.activeSpecies();
    const basicSpecieMixture& compCont = thermoCont.composition();

    // Set global species indices

    labelList indices(2);

    forAll(species, j)
    {
        indices[j] = activeSpecies[species[j]];
    }

    const scalarField sumZ(min(thermo_.sumZ()->field(),1.0));

    PtrList<scalarField> gamma(2);

    // Initialize gamma as 1

    forAll(species, j)
    {
        gamma.set
        (
            j,
            new scalarField(sumZ.size(),1.0)
        );
    }

    scalarList W(2);

    forAll(species, j)
    {
        W[j] = compCont.W(indices[j]);
    }

    forAll(sumZ, celli)
    {
        if (sumZ[celli] > VSMALL)
        {
            scalarList Z(2);

            forAll(species, j)
            {
                Z[j] = thermo_.Z()[indices[j]][celli];
            }

            const scalarList z(Z/sumZ[celli]);
            const scalarList w(z/W/sum(z/W));

            const label jA(thermo_.activeSpecies()[firstSpecieName_]);
            const label jB(jA == 0 ? 1 : 0);

            gamma[jA][celli] =
                Foam::exp
                (
                    C1_/(1.0 + C1_/C2_*(w[jA]/max(1.0-w[jA],VSMALL)))
                );

            gamma[jB][celli] =
                Foam::exp
                (
                    C2_/C1_
                  * (
                        2.0*Foam::sqrt(C1_*Foam::log(gamma[jA][celli]))
                      + Foam::log(gamma[jA][celli])
                      + C1_
                    )
                );
        }
    }

    return gamma;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
