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

#include "dispersedDiffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dispersedDiffusionModel> dispersedDiffusionModel::New
(
    aerosolModel& aerosol,
    const dictionary& dict
)
{
    const word modelType(dict.lookup("type"));

    Info<< "    Selecting dispersed diffusion model " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown dispersedDiffusion model type " << modelType
            << nl << nl << "Valid dispersedDiffusionModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return
        autoPtr<dispersedDiffusionModel>
        (
            cstrIter()(aerosol, dict)
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
