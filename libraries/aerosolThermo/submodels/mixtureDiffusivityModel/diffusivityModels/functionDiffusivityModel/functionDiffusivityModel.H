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

#ifndef functionDiffusivityModel_H
#define functionDiffusivityModel_H

#include "diffusivityModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

/*---------------------------------------------------------------------------*\
                   Class functionDiffusivityModel Declaration
\*---------------------------------------------------------------------------*/

class functionDiffusivityModel
:
    public diffusivityModel
{
    // Private data

        autoPtr<Function1<scalar>> function_;

    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);

        //- Disallow default bitwise assignment
        void operator=(const functionDiffusivityModel&);


public:

    // Runtime type information
    TypeName("function");


    // Constructors

        //- Construct from entry name and dictionary
        functionDiffusivityModel
        (
            const word& entryName,
            const dictionary& dict,
            aerosolThermo& thermo,
            const label j,
            const label k
        );

        //- Copy constructor
        functionDiffusivityModel(const functionDiffusivityModel& model);

        //- Construct and return a clone
        virtual tmp<diffusivityModel> clone() const
        {
            return tmp<diffusivityModel>(new functionDiffusivityModel(*this));
        }


    //- Destructor
    virtual ~functionDiffusivityModel();


    // Member Functions

        //- Return the binary diffusivity field
        virtual tmp<scalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
