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

\file componentMixingModel.H
\brief Base class for the phase mixing model

The componentMixingModel object determines how a property of both the continuous
and the dispersed phase is mixed into a mixture property. The mixing can be
done for a volField, cell list or patch.

*/

#ifndef componentMixingModel_H
#define componentMixingModel_H

#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward decleration
class rhoAerosolPhaseThermo;

/*---------------------------------------------------------------------------*\
                       Class componentMixingModel Declaration
\*---------------------------------------------------------------------------*/

class componentMixingModel
{
public:

    // Public data types

        //- Thermo field type enumeration
        enum fieldType
        {
            viscosity,
            conductivity
        };

        //- Thermo field type names
        static const Enum<fieldType> fieldTypeNames_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        componentMixingModel(const componentMixingModel&);

        //- Disallow default bitwise assignment
        void operator=(const componentMixingModel&);


protected:

    // Protected Data

        //- Reference to the phase thermo object
        const rhoAerosolPhaseThermo& thermo_;

        //- Type of field to mix
        fieldType field_;


public:

    //- Runtime type information
    TypeName("componentMixingModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             componentMixingModel,
             dictionary,
             (
                const rhoAerosolPhaseThermo& thermo,
                const fieldType& mixingModel
             ),
             (thermo, mixingModel)
         );

    // Constructors

        //- Construct from type name, dictionary and thermo model
        componentMixingModel
        (
            const word& modelType,
            const rhoAerosolPhaseThermo& thermo,
            const fieldType& field
        );


    // Selectors

        //- Return a reference to the selected mixing model
        static autoPtr<componentMixingModel> New
        (
            const word& modelType,
            const rhoAerosolPhaseThermo& thermo,
            const fieldType& field
        );


    //- Destructor
    virtual ~componentMixingModel();


    // Member Functions

        // Access

            //- Reference to the phase thermo object
            inline const rhoAerosolPhaseThermo& thermo()
            {
                return thermo_;
            }

            //- Return the mixing model type
            inline const fieldType& field()
            {
                return field_;
            }

        // Mixing

            //- Set the mixing at the given field
            virtual void mix(volScalarField& phi) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
