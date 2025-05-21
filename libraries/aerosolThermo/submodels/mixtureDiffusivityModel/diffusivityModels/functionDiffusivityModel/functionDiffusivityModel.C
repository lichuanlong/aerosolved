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

#include "functionDiffusivityModel.H"
#include "aerosolThermo.H"
#include "makeDiffusivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{
    makeDiffusivityModel(functionDiffusivityModel, diffusivityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::diffusivityModels::functionDiffusivityModel::read
(
    const dictionary& coeffs
)
{}


Foam::diffusivityModels::functionDiffusivityModel::functionDiffusivityModel
(
    const word& entryName,
    const dictionary& dict,
    aerosolThermo& thermo,
    const label j,
    const label k
)
:
    diffusivityModel(entryName, thermo, j, k),
    function_()
{
    Info<< dict << endl;
    Info<< entryName << endl;
    function_ =
        Function1<scalar>::New(entryName, dict);
    Info<< "nice2" << endl;
    read(dict);
}


Foam::diffusivityModels::functionDiffusivityModel::functionDiffusivityModel
(
    const functionDiffusivityModel& model
)
:
    diffusivityModel(model),
    function_(model.function_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusivityModels::functionDiffusivityModel::~functionDiffusivityModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::diffusivityModels::functionDiffusivityModel::D() const
{
    return function_->value(thermo_.thermoCont().T().field());
}

// ************************************************************************* //
