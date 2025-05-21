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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<activityCoeffModel> activityCoeffModel::New(aerosolThermo& thermo)
{
    word modelType;

    if (thermo.found(activityCoeffModel::dictName))
    {
        modelType = word
        (
            thermo.subDict(activityCoeffModel::dictName).lookup("type")
        );
    }
    else
    {
        modelType = "none";
    }

    Info<< "    Selecting activity coefficient model " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown activity coefficient model type " << modelType
            << nl << nl << "Valid activity coefficient model types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return
        autoPtr<activityCoeffModel>
        (
            cstrIter()(thermo)
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
