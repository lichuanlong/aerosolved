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

#ifndef WilkeWassiljewaComponentMixing_H
#define WilkeWassiljewaComponentMixing_H

#include "componentMixingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WilkeWassiljewaComponentMixing Declaration
\*---------------------------------------------------------------------------*/

class WilkeWassiljewaComponentMixing
:
    public componentMixingModel
{
    // Private member functions

        //- Disallow WilkeWassiljewa bitwise copy construct
        WilkeWassiljewaComponentMixing(const WilkeWassiljewaComponentMixing&);

        //- Disallow WilkeWassiljewa bitwise assignment
        void operator=(const WilkeWassiljewaComponentMixing&);


public:

    //- Runtime type information
    TypeName("WilkeWassiljewa");


    // Constructors

        //- Construct from thermo model
        WilkeWassiljewaComponentMixing
        (
            const rhoAerosolPhaseThermo& thermo,
            const fieldType& field
        );


    //- Destructor
    virtual ~WilkeWassiljewaComponentMixing();


    // Member Functions

        //- Set the mixing at the given field
        virtual void mix(volScalarField& phi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
