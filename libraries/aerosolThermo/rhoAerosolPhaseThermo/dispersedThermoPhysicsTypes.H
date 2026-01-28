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

#ifndef dispersedThermoPhysicsTypes_H
#define dispersedThermoPhysicsTypes_H

#include "specie.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"
#include "rhoConst.H"
#include "eConstThermo.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "icoFunction.H"

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
                perfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    > constDispHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    > constDispEHThermoPhysicsAeroSolved;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    > dispHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                icoFunction<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleFuncDispHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                icoFunction<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleFuncDispEHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                icoPolynomial<specie, 8>
            >,
            sensibleEnthalpy
        >
    > constIncompressiblePoly8DispHThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                icoPolynomial<specie, 8>
            >,
            sensibleEnthalpy
        >
    > constIncompressiblePoly8DispEHThermoPhysicsAeroSolved;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    > adiabaticDispHThermoPhysicsAeroSolved;

    // thermo physics types based on sensibleInternalEnergy

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    > constDispEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    > constDispHEThermoPhysicsAeroSolved;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    > dispEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                icoFunction<specie>
            >,
            sensibleInternalEnergy
        >
    > constIncompressibleFuncDispEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                icoFunction<specie>
            >,
            sensibleInternalEnergy
        >
    > constIncompressibleFuncDispHEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                icoPolynomial<specie, 8>
            >,
            sensibleInternalEnergy
        >
    > constIncompressiblePoly8DispEThermoPhysicsAeroSolved;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                icoPolynomial<specie, 8>
            >,
            sensibleInternalEnergy
        >
    > constIncompressiblePoly8DispHEThermoPhysicsAeroSolved;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    > adiabaticDispEThermoPhysicsAeroSolved;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
