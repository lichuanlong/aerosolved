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

if (adjustTimeStep)
{
    // By default, take maxCo as the max CFL number for the aerosol source term,
    // drift and corrective drift. Can be overwritten by specifying
    // maxAerosolCo, maxDriftCo and maxCorrCo in the controlDict

    const scalar maxAerosolCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxAerosolCo", maxCo);

    const scalar maxDriftCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxDriftCo", maxCo);

    const scalar maxCorrCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCorrCo", maxCo);

    const scalar maxConvDeltaTFact = maxCo/(CoNum + SMALL);
    const scalar maxAerosolDeltaTFact = maxAerosolCo/(aerosolCoNum + SMALL);
    const scalar maxDriftDeltaTFact = maxDriftCo/(driftCoNum + SMALL);
    const scalar maxCorrDeltaTFact = maxCorrCo/(corrCoNum + SMALL);

    scalar maxDeltaTFact = 1e99;

    // Take the maxDeltaTFact as the minimum of all constraints

    maxDeltaTFact = min(maxDeltaTFact,maxConvDeltaTFact);
    maxDeltaTFact = min(maxDeltaTFact,maxAerosolDeltaTFact);
    maxDeltaTFact = min(maxDeltaTFact,maxDriftDeltaTFact);
    maxDeltaTFact = min(maxDeltaTFact,maxCorrDeltaTFact);

    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}
