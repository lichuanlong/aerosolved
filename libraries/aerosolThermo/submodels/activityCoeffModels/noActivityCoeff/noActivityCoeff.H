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

#ifndef noActivityCoeff_H
#define noActivityCoeff_H

#include "activityCoeffModel.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noActivityCoeff Declaration
\*---------------------------------------------------------------------------*/

class noActivityCoeff
:
    public activityCoeffModel
{

    // Private member functions

        //- Disallow default bitwise copy construct
        noActivityCoeff(const noActivityCoeff&);

        //- Disallow default bitwise assignment
        void operator=(const noActivityCoeff&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from thermo model
        noActivityCoeff(aerosolThermo& thermo);


    //- Destructor
    virtual ~noActivityCoeff();


    // Member Functions

        // Evolution

            //- Compute the activity coefficients for a list of species
            virtual PtrList<scalarField> gamma
            (
                const speciesTable& species
            ) const;

            //- Compute the activity coefficient for a list of species andgiven
            //  pressure, temperature and dispersed molar composition fields
            virtual PtrList<scalarField> gamma
            (
                const speciesTable& species,
                const scalarField& p,
                const scalarField& T,
                const PtrList<scalarField>& x
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
