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

\file componentMixing.H
\brief Class which contains the three componentMixingModels

The componentMixingModels are for viscosity, conductivity and heat capacity

*/

#ifndef componentMixing_H
#define componentMixing_H

#include "autoPtr.H"
#include "componentMixingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class componentMixing Declaration
\*---------------------------------------------------------------------------*/

class componentMixing
{
private:

    // Private Data

        //- Reference to the aerosol phase thermo object
        const rhoAerosolPhaseThermo& thermo_;

        //- Viscosity mixing object
        autoPtr<componentMixingModel> viscosity_;

        //- Conductivity mixing object
        autoPtr<componentMixingModel> conductivity_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        componentMixing(const componentMixing&);

        //- Disallow default bitwise assignment
        void operator=(const componentMixing&);


public:

    // Constructors

        //- Construct from type name, dictionary and thermo model
        componentMixing
        (
            const rhoAerosolPhaseThermo& thermo
        );


    //- Destructor
    virtual ~componentMixing();


    // Member Functions

        // Access

            //- Reference to the phase thermo object
            inline const rhoAerosolPhaseThermo& thermo()
            {
                return thermo_;
            }

        // Mixing

            //- Mixing of viscosity
            inline componentMixingModel& viscosity()
            {
                return viscosity_();
            }

            //- Mixing of conductivity
            inline componentMixingModel& conductivity()
            {
                return conductivity_();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
