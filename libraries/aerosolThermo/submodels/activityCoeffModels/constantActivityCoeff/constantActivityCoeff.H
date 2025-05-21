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

#ifndef constantActivityCoeff_H
#define constantActivityCoeff_H

#include "activityCoeffModel.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constantActivityCoeff Declaration
\*---------------------------------------------------------------------------*/

class constantActivityCoeff
:
    public activityCoeffModel
{
    // Private Data

        //- List of constant coefficients read from dictionary
        scalarList coeffs_;


    // Private member functions

        //- Disallow default bitwise copy construct
        constantActivityCoeff(const constantActivityCoeff&);

        //- Disallow default bitwise assignment
        void operator=(const constantActivityCoeff&);


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from thermo model
        constantActivityCoeff(aerosolThermo& thermo);


    //- Destructor
    virtual ~constantActivityCoeff();


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
