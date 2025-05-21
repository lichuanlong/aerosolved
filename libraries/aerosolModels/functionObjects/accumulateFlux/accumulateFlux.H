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

\file accumulateFlux.H

\brief functionObject to time-integrate a given flux in order to compute the
accumulated flux. The flux is scaled by the surface area of each cell, to
compute a deposition density in units of [flux]*s/m2 where [flux] is the unit of
the specified flux field.

*/

#ifndef functionObjects_accumulateFlux_H
#define functionObjects_accumulateFlux_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class accumulateFlux Declaration
\*---------------------------------------------------------------------------*/

class accumulateFlux
:
    public fieldExpression
{
    // Private Data

        //- Accumulation field pointer
        autoPtr<volScalarField> accumulatedPtr_;

    // Private Member Functions

        //- Update the accumulated flux
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("accumulateFlux");


    // Constructors

        //- Construct from Time and dictionary
        accumulateFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Called at each ++ or += of the time-loop. Overwrites default
        //  behavior based on the executeControl parameter
        virtual bool execute();

    //- Destructor
    virtual ~accumulateFlux();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
