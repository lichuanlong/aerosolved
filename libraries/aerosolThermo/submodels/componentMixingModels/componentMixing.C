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

#include "componentMixing.H"
#include "rhoAerosolPhaseThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define viscMixModel componentMixingModel::fieldType::viscosity
#define condMixModel componentMixingModel::fieldType::conductivity

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

componentMixing::componentMixing
(
    const rhoAerosolPhaseThermo& thermo
)
:
    thermo_(thermo),
    viscosity_(),
    conductivity_()
{
    word viscMixing("default");
    word condMixing("default");

    const rhoThermo& rthermo = thermo;

    if (rthermo.found("componentMixing"))
    {
        const dictionary& compMix = rthermo.subDict("componentMixing");

        viscMixing =
            compMix.lookupOrDefault<word>("viscosity", viscMixing);

        condMixing =
            compMix.lookupOrDefault<word>("conductivity", condMixing);
    }

    viscosity_ =
        componentMixingModel::New(viscMixing, thermo, viscMixModel);

    conductivity_ =
        componentMixingModel::New(condMixing, thermo, condMixModel);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

componentMixing::~componentMixing()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
