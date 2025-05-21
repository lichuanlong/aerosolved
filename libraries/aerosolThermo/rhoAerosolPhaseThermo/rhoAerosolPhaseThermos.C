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

#include "rhoAerosolPhaseThermo.H"
#include "heAerosolRhoThermo.H"
#include "heRhoThermo.H"

#include "makeAerosolThermo.H"
#include "addToRunTimeSelectionTable.H"

#include "dispersedThermoPhysicsTypes.H"
#include "continuousThermoPhysicsTypes.H"
#include "aerosolPhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Continuous models, sensibleEnthalpy

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constGasHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constGasEHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    gasHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressibleGasHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressibleGasEHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    incompressibleGasHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    icoPoly8HThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    icoPoly8PerfectGasHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    icoPoly8RhoConstHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constAdiabaticFluidHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constAdiabaticFluidEHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constEHThermoPhysicsAeroSolved
);

// Continuous models, sensibleInternalEnergy

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constGasEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constGasHEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    gasEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressibleGasEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressibleGasHEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    incompressibleGasEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constAdiabaticFluidEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constAdiabaticFluidHEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constHEThermoPhysicsAeroSolved
);

// Dispersed models, sensibleEnthalpy

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constDispHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constDispEHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    dispHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressibleFuncDispHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressiblePoly8DispHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressiblePoly8DispEHThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    adiabaticDispHThermoPhysicsAeroSolved
);

// Dispersed models, sensibleInternalEnergy

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constDispEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constDispHEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    dispEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressiblePoly8DispEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressiblePoly8DispHEThermoPhysicsAeroSolved
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    adiabaticDispEThermoPhysicsAeroSolved
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
