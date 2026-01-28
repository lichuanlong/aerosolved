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

#include "HirschfelderCurtiss.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(HirschfelderCurtiss, 0);
addToRunTimeSelectionTable(continuousDiffusionModel, HirschfelderCurtiss, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HirschfelderCurtiss::HirschfelderCurtiss
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    continuousDiffusionModel(type(), aerosol, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HirschfelderCurtiss::~HirschfelderCurtiss()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> HirschfelderCurtiss::D(const label j)
{
    aerosolThermo& thermo = aerosol_.thermo();

    tmp<volScalarField> tD
    (
        new volScalarField
        (
            IOobject
            (
                word(thermo.Y()[j].name() + ":D"),
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aerosol_.mesh(),
            dimensionedScalar("D", sqr(dimLength)/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& D = tD.ref();

    D.field() = thermo.diffusivity().Deff(j);
    D.correctBoundaryConditions();

    return D*(1.0 - thermo.Y()[j]/thermo.sumY());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
