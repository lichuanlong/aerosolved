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

#ifndef noDispersedInertialDrift_H
#define noDispersedInertialDrift_H

#include "dispersedInertialDriftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDispersedInertialDrift Declaration
\*---------------------------------------------------------------------------*/

class noDispersedInertialDrift
:
    public dispersedInertialDriftModel
{
private:

    // Private Member Data

        // Dummy field to return
        const volVectorField dummy_;

    // Private member functions

        //- Disallow default bitwise copy construct
        noDispersedInertialDrift(const noDispersedInertialDrift&);

        //- Disallow default bitwise assignment
        void operator=(const noDispersedInertialDrift&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from aerosol model
        noDispersedInertialDrift(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~noDispersedInertialDrift();


    // Member Functions

        //- Compute and return the relative velocity given a droplet size field
        //  and size name
        virtual const volVectorField& V
        (
            const volScalarField& d,
            const word sizeName
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
