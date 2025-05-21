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

\file boundaryFluxDensity.H
\brief functionObject that reads a flux field and outputs the boundary flux
density as a volField

*/

#ifndef functionObjects_boundaryFluxDensity_H
#define functionObjects_boundaryFluxDensity_H

#include "fvMeshFunctionObject.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class boundaryFluxDensity Declaration
\*---------------------------------------------------------------------------*/

class boundaryFluxDensity
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of the flux field
        word phiName_;

        //- Name of the result field
        word resultName_;


    // Private Member Functions

        //- Calculate the diameter field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("boundaryFluxDensity");


    // Constructors

        //- Construct from Time and dictionary
        boundaryFluxDensity
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~boundaryFluxDensity();

    //- Member Functions
        //- Calculate the result field
        virtual bool execute();

        //- Write the result field
        virtual bool write();

        //- Clear the result field from the objectRegistry
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
