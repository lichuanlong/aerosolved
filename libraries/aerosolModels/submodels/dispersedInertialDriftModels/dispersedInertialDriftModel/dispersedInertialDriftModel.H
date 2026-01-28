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

\file dispersedInertialDriftModel.H
\brief Base class of the dispersed inertial drift model

*/

#ifndef dispersedInertialDriftModel_H
#define dispersedInertialDriftModel_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dispersedInertialDriftModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedInertialDriftModel
:
    public aerosolSubModelBase
{
private:

    // Private Data

        //- Drift velocity base field, needed for boundary conditions
        autoPtr<volVectorField> VPtr_;

        //- Drift velocity cache
        HashPtrTable<volVectorField> cache_;

        //- Max relative velocity (defaults to 1e99)
        scalar VMax_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dispersedInertialDriftModel(const dispersedInertialDriftModel&);

        //- Disallow default bitwise assignment
        void operator=(const dispersedInertialDriftModel&);


protected:

    // Protected Data

        //- Gravitational accelleration
        uniformDimensionedVectorField g_;

    // Protected Member Functions

        //- Return a cached field. If it does not exist, create it.
        volVectorField& getCachedField(const word sizeName);

        //- Function to read the base drift velocity field. To be called by
        //  derived types.
        void readBaseField();

        //- Function to read the base drift velocity field if present, and
        //  otherwise create a dummy one. To be called by derived types.
        void readBaseFieldIfPresent();

public:

    //- Runtime type information
    TypeName("dispersedInertialDriftModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             dispersedInertialDriftModel,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict
             ),
             (aerosol, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        dispersedInertialDriftModel
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected dispersed inertial drift model
        static autoPtr<dispersedInertialDriftModel> New
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dispersedInertialDriftModel();


    // Member Functions

        //- Return the velocity field name given a size name
        virtual word velocityFieldName(const word sizeName) const;

        //- Return the drift velocity cache
        inline HashPtrTable<volVectorField>& cache()
        {
            return cache_;
        }

        //- Compute the particle Reynolds number
        tmp<volScalarField> Re
        (
            const volScalarField& d,
            const volVectorField& V
        ) const;

        //- Function to limit the velocity
        virtual void limit(volVectorField& V) const;

        //- Compute and return the relative velocity given a droplet size field
        //  and size name
        virtual const volVectorField& V
        (
            const volScalarField& d,
            const word sizeName
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
