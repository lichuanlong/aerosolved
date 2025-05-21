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

\file driftFluxModel.H
\brief Container class of the drift flux modeling components

This class contains the continuous phase diffusion, dispersed phase diffusion
and dispersed inertial drift objects, as well as fluxes and multivariate
convection schemes. The class is mostly a container, creating some layering and
structure in the code. Its objects are manipulated by the aerosol model.

*/

#ifndef driftFluxModel_H
#define driftFluxModel_H

#include "aerosolSubModelBase.H"
#include "continuousDiffusionModel.H"
#include "dispersedDiffusionModel.H"
#include "dispersedInertialDriftModel.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class driftFluxModel Declaration
\*---------------------------------------------------------------------------*/

class driftFluxModel
:
    public aerosolSubModelBase
{
private:

    // Private Member Data

        //- Continuous diffusion model
        autoPtr<continuousDiffusionModel> contDiff_;

        //- Dispersed diffusion model
        autoPtr<dispersedDiffusionModel> dispDiff_;

        //- Dispersed inertial drift model
        autoPtr<dispersedInertialDriftModel> dispInertialDrift_;

        //- Multivariate convection scheme for phi
        tmp<fv::convectionScheme<scalar>> mvPhi_;

        //- Multivariate convection scheme for phiCorr
        tmp<fv::convectionScheme<scalar>> mvPhiCorr_;

        //- Multivariate convection scheme for phiInertial
        tmp<fv::convectionScheme<scalar>> mvPhiInertial_;

        //- Table of fields for multivariate convection
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;

        //- Corrective mass fraction drift flux [kg/s]
        surfaceScalarField phiCorr_;

        //- Corrective mass fraction drift flux for diffusion [kg/s]
        surfaceScalarField phiCorrDiff_;

        //- Dispersed inertial drift mass fraction drift flux [kg/s]
        surfaceScalarField phiInertial_;

        //- Dispersed mass fraction diffusivity [m^2/s]
        volScalarField DDisp_;

        //- Continuous mass fraction diffusivities [m^2/s]
        PtrList<volScalarField> DCont_;

        //- Relative phase velocity [m/s]
        volVectorField Ur_;

        //- Continuous phase velocity [m/s]
        volVectorField Uc_;

        //- Dispersed phase velocity [m/s]
        volVectorField Ud_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        driftFluxModel(const driftFluxModel&);

        //- Disallow default bitwise assignment
        void operator=(const driftFluxModel&);


public:

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        driftFluxModel
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~driftFluxModel();


    // Access

        //- Access to the continuous diffusion model
        virtual inline continuousDiffusionModel& contDiff()
        {
            return contDiff_();
        }

        //- Const access to the continuous diffusion model
        virtual inline const continuousDiffusionModel& contDiff() const
        {
            return contDiff_();
        }

        //- Access to the dispersed diffusion model
        virtual inline dispersedDiffusionModel& dispDiff()
        {
            return dispDiff_();
        }

        //- Const access to the dispersed diffusion model
        virtual inline const dispersedDiffusionModel& dispDiff() const
        {
            return dispDiff_();
        }

        //- Access to the dispersed inertial drift model
        virtual inline dispersedInertialDriftModel& dispInertialDrift()
        {
            return dispInertialDrift_();
        }

        //- Const access to the dispersed inertial drift model
        virtual inline const dispersedInertialDriftModel& dispInertialDrift() const
        {
            return dispInertialDrift_();
        }

        //- Return reference to multivariate convectionScheme for phi
        inline const fv::convectionScheme<scalar>& mvPhi() const;

        //- Return reference to multivariate convectionScheme for phiCorr
        inline const fv::convectionScheme<scalar>& mvPhiCorr() const;

        //- Return reference to multivariate convectionScheme for phiInertial
        inline const fv::convectionScheme<scalar>& mvPhiInertial() const;

        //- Return reference to the able of fields for multivariate convection
        inline multivariateSurfaceInterpolationScheme<scalar>::fieldTable&
        fields();

        //- Return reference to the able of fields for multivariate convection
        inline const multivariateSurfaceInterpolationScheme<scalar>::fieldTable&
        fields() const;

        //- Return the relative phase velocity
        inline volVectorField& Ur();

        //- Return the relative phase velocity
        inline const volVectorField& Ur() const;

        //- Return the continuous phase velocity
        inline const volVectorField& Uc() const;

        //- Return the dispersed phase velocity
        inline const volVectorField& Ud() const;

        //- Return non-const access to the corrective mass fraction drift flux
        inline surfaceScalarField& phiCorr();

        //- Return access to the corrective mass fraction drift flux
        inline const surfaceScalarField& phiCorr() const;

        //- Return non-const access to the corrective mass fraction drift flux
        //  for diffusion
        inline surfaceScalarField& phiCorrDiff();

        //- Return access to the corrective mass fraction drift flux for
        //  diffusion
        inline const surfaceScalarField& phiCorrDiff() const;

        //- Return access to the dispersed inertial mass fraction drift flux
        inline surfaceScalarField& phiInertial();

        //- Return access to the dispersed inertial mass fraction drift flux
        inline const surfaceScalarField& phiInertial() const;

        //- Return access to the dispersed mass fraction diffusivity
        inline volScalarField& DDisp();

        //- Return access to the dispersed mass fraction diffusivity
        inline const volScalarField& DDisp() const;

        //- Return access to the continuous mass fraction diffusivities
        inline PtrList<volScalarField>& DCont();

        //- Return access to the continuous mass fraction diffusivities
        inline const PtrList<volScalarField>& DCont() const;

    // Member Functions

        //- Correct the phase velocities
        virtual void correctPhaseVelocities();

        //- Correct the drift flux models
        virtual void correct();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "driftFluxModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
