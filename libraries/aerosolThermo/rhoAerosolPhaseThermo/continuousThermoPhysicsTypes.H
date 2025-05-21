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

#ifndef continuousThermoPhysicsTypes_H
#define continuousThermoPhysicsTypes_H

#include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"
#include "rhoConst.H"
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include "sutherlandTransport.H"
#include "constTransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // thermo physics types based on sensibleEnthalpy

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constGasHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constGasEHThermoPhysicsAeroSolved;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > gasHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleGasHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleGasEHThermoPhysicsAeroSolved;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > incompressibleGasHThermoPhysicsAeroSolved;

    typedef
    polynomialTransport
    <
        species::thermo
        <
            hPolynomialThermo
            <
                icoPolynomial<specie, 8>,
                8
            >,
            sensibleEnthalpy
        >,
        8
    > icoPoly8HThermoPhysicsAeroSolved;

    typedef
    polynomialTransport
    <
        species::thermo
        <
            hPolynomialThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >,
        8
    > icoPoly8PerfectGasHThermoPhysicsAeroSolved;

    typedef
    polynomialTransport
    <
        species::thermo
        <
            hPolynomialThermo
            <
                rhoConst<specie>
            >,
            sensibleEnthalpy
        >,
        8
    > icoPoly8RhoConstHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    >
    constAdiabaticFluidHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    >
    constAdiabaticFluidEHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                rhoConst<specie>
            >,
            sensibleEnthalpy
        >
    >
    constHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                rhoConst<specie>
            >,
            sensibleEnthalpy
        >
    >
    constEHThermoPhysicsAeroSolved;

    // thermo physics types based on sensibleInternalEnergy

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > constGasEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > constGasHEThermoPhysicsAeroSolved;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > gasEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > constIncompressibleGasEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > constIncompressibleGasHEThermoPhysicsAeroSolved;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > incompressibleGasEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constAdiabaticFluidEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constAdiabaticFluidHEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                rhoConst<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                rhoConst<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constHEThermoPhysicsAeroSolved;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
