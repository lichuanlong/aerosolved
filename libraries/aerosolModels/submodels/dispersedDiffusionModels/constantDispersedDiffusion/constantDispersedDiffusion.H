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

#ifndef constantDispersedDiffusion_H
#define constantDispersedDiffusion_H

#include "dispersedDiffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constantDispersedDiffusion Declaration
\*---------------------------------------------------------------------------*/

class constantDispersedDiffusion
:
    public dispersedDiffusionModel
{
private:

    // Private Data

        //- Constant diffusivity value
        scalar diffusivity_;

    // Private member functions

        //- Disallow default bitwise copy construct
        constantDispersedDiffusion(const constantDispersedDiffusion&);

        //- Disallow default bitwise assignment
        void operator=(const constantDispersedDiffusion&);

        //- Compute the diffusivity given droplet size, pressure and
        //  temperature
        scalar diffusivity
        (
            const scalar& d,
            const scalar& p,
            const scalar& T
        ) const;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from aerosol model
        constantDispersedDiffusion(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~constantDispersedDiffusion();


    // Member Functions

        //- Compute the diffusivity given a droplet size and cell
        scalar D(const scalar& d, const label& celli) const;

        //- Compute the diffusivity given a droplet size field and patch
        scalarField D(const scalarField& d, const label& patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
