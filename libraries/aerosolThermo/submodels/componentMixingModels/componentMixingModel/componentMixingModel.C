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

#include "componentMixingModel.H"
#include "rhoAerosolPhaseThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(componentMixingModel, 0);
defineRunTimeSelectionTable(componentMixingModel, dictionary);

const Enum
<
    componentMixingModel::fieldType
>
componentMixingModel::fieldTypeNames_
{
    { fieldType::viscosity, "viscosity" },
    { fieldType::conductivity, "conductivity" },
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

componentMixingModel::componentMixingModel
(
    const word& modelType,
    const rhoAerosolPhaseThermo& thermo,
    const fieldType& field
)
:
    thermo_(thermo),
    field_(field)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

componentMixingModel::~componentMixingModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
