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

#ifndef defaultComponentMixing_H
#define defaultComponentMixing_H

#include "componentMixingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class defaultComponentMixing Declaration
\*---------------------------------------------------------------------------*/

class defaultComponentMixing
:
    public componentMixingModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        defaultComponentMixing(const defaultComponentMixing&);

        //- Disallow default bitwise assignment
        void operator=(const defaultComponentMixing&);


public:

    //- Runtime type information
    TypeName("default");


    // Constructors

        //- Construct from thermo model
        defaultComponentMixing
        (
            const rhoAerosolPhaseThermo& thermo,
            const fieldType& field
        );


    //- Destructor
    virtual ~defaultComponentMixing();


    // Member Functions

        //- Set the mixing at the given field
        virtual void mix(volScalarField& phi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
