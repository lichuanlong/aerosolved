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

\file sectionalPecletNumbers.H
\brief functionObject to compute the sectional Peclet numbers at faces

*/

#ifndef functionObjects_sectionalPecletNumbers_H
#define functionObjects_sectionalPecletNumbers_H

#include "PecletNumber.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class sectionalPecletNumbers Declaration
\*---------------------------------------------------------------------------*/

class sectionalPecletNumbers
:
    public PecletNumber
{
    // Private Member Functions

        //- Calculate the sectional Peclet number fields and return true if
        //  successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("sectionalPecletNumbers");


    // Constructors

        //- Construct from Time and dictionary
        sectionalPecletNumbers
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sectionalPecletNumbers();

    //- Member Functions

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
