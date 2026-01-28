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

/**

\file activityCoeffModel.H
\brief Base class for the activity coefficient model

*/

#ifndef activityCoeffModel_H
#define activityCoeffModel_H

#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "PtrList.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class activityCoeffModel Declaration
\*---------------------------------------------------------------------------*/

//Forward declaration
class aerosolThermo;

class activityCoeffModel
{
protected:

    // Protected Member Data

        //- Reference to the thermo model
        aerosolThermo& thermo_;

        //- Reference to the dictionary
        const dictionary& dict_;


    // Protected Member Functions

        //- Unity gamma fields of given size
        PtrList<scalarField> unityGamma
        (
            const label nSpecies,
            const label nCells
        ) const;

        //- Sum the given list of fields
        tmp<scalarField> sumList(const PtrList<scalarField>& x) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        activityCoeffModel(const activityCoeffModel&);

        //- Disallow default bitwise assignment
        void operator=(const activityCoeffModel&);


public:

    //- Runtime type information
    TypeName("activityCoeffModel");

    //- Name of the activity coefficient dictionary
    static const word dictName;


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             activityCoeffModel,
             dictionary,
             (aerosolThermo& thermo),
             (thermo)
         );

    // Constructors

        //- Construct from thermo
        activityCoeffModel(aerosolThermo& thermo);


    // Selectors

        //- Return a reference to the selected activity coeff model
        static autoPtr<activityCoeffModel> New(aerosolThermo& thermo);


    //- Destructor
    virtual ~activityCoeffModel();


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
            ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
