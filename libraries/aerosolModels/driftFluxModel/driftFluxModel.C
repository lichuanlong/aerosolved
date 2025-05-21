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

#include "driftFluxModel.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

driftFluxModel::driftFluxModel
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    aerosolSubModelBase
    (
        aerosol,
        dict.subDict("driftFluxModel"),
        "driftFluxModel",
        "driftFluxModel"
    ),
    contDiff_(),
    dispDiff_(),
    dispInertialDrift_(),
    phiCorr_
    (
        IOobject
        (
            "phiCorr",
            aerosol_.mesh().time().timeName(),
            aerosol_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux
        (
            volVectorField
            (
                IOobject
                (
                    "dummy",
                    aerosol_.mesh().time().timeName(),
                    aerosol_.mesh()
                ),
                aerosol_.mesh(),
                dimensionedVector
                (
                    "dummy",
                    dimVelocity*dimDensity,
                    vector::zero
                )
            )
        )
    ),
    phiCorrDiff_
    (
        IOobject
        (
            "phiCorrDiff",
            aerosol_.mesh().time().timeName(),
            aerosol_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux
        (
            volVectorField
            (
                IOobject
                (
                    "dummy",
                    aerosol_.mesh().time().timeName(),
                    aerosol_.mesh()
                ),
                aerosol_.mesh(),
                dimensionedVector
                (
                    "dummy",
                    dimVelocity*dimDensity,
                    vector::zero
                )
            )
        )
    ),
    phiInertial_
    (
        IOobject
        (
            "phiInertial",
            aerosol_.mesh().time().timeName(),
            aerosol_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux
        (
            volVectorField
            (
                IOobject
                (
                    "dummy",
                    aerosol_.mesh().time().timeName(),
                    aerosol_.mesh()
                ),
                aerosol_.mesh(),
                dimensionedVector
                (
                    "dummy",
                    dimVelocity*dimDensity,
                    vector::zero
                )
            )
        )
    ),
    DDisp_
    (
        IOobject
        (
            "DDisp",
            aerosol_.mesh().time().timeName(),
            aerosol_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aerosol_.mesh(),
        dimensionedScalar("D", dimArea/dimTime, 0.0)
    ),
    DCont_(aerosol_.thermo().contSpecies().size()),
    Ur_
    (
        IOobject
        (
            "Ur",
            aerosol_.mesh().time().timeName(),
            aerosol_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        aerosol_.mesh(),
        dimensionedVector("Ur", dimVelocity, vector::zero)
    ),
    Uc_
    (
        IOobject
        (
            "Uc",
            aerosol_.mesh().time().timeName(),
            aerosol_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        aerosol_.mesh(),
        dimensionedVector("Uc", dimVelocity, vector::zero)
    ),
    Ud_
    (
        IOobject
        (
            "Ud",
            aerosol_.mesh().time().timeName(),
            aerosol_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        aerosol_.mesh(),
        dimensionedVector("Ud", dimVelocity, vector::zero)
    )
{
    // Create the continuous phase diffusion model from the continuousDiffusion
    // dictionary. Add legacy support for the diffusion dictionary

    contDiff_ =
        continuousDiffusionModel::New
        (
            aerosol,
            this->dict().found("diffusion")
        && !this->dict().found("continuousDiffusion")
          ? this->dict().subDict("diffusion")
          : this->dict().subDict("continuousDiffusion")
        );

    // Create the dispersed phase diffusion model from the dispersedDiffusion
    // dictionary. Add legacy support for the Brownian dictionary

    dispDiff_ =
        dispersedDiffusionModel::New
        (
            aerosol,
            this->dict().found("Brownian")
        && !this->dict().found("dispersedDiffusion")
          ? this->dict().subDict("Brownian")
          : this->dict().subDict("dispersedDiffusion")
        );

    // Create the dispersed phase inertial drift model from the
    // dispersedInertialDrift dictionary. Add legacy support for the inertial
    // dictionary

    dispInertialDrift_ =
        dispersedInertialDriftModel::New
        (
            aerosol,
            this->dict().found("inertial")
        && !this->dict().found("dispersedInertialDrift")
          ? this->dict().subDict("inertial")
          : this->dict().subDict("dispersedInertialDrift")
        );

    forAll(aerosol_.thermo().contSpecies(), j)
    {
        DCont_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("DCont", aerosol_.thermo().contSpecies()[j]),
                    aerosol_.mesh().time().timeName(),
                    aerosol_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aerosol_.mesh(),
                dimensionedScalar("DCont", dimArea/dimTime, 0.0)
            )
        );
    }

    forAll(aerosol_.thermo().contSpecies(), j)
    {
        fields_.add(aerosol_.thermo().Y()[j]);
    }

    forAll(aerosol_.thermo().dispSpecies(), j)
    {
        fields_.add(aerosol_.thermo().Z()[j]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

driftFluxModel::~driftFluxModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void driftFluxModel::correct()
{
    // Set the multivariate convection schemes

    mvPhi_ =
        fv::convectionScheme<Foam::scalar>::New
        (
            aerosol_.mesh(),
            fields_,
            aerosol_.phi(),
            aerosol_.mesh().divScheme("div(mvConv)")
        );

    mvPhiCorr_ =
        fv::convectionScheme<Foam::scalar>::New
        (
            aerosol_.mesh(),
            fields_,
            phiCorr_,
            aerosol_.mesh().divScheme("div(mvConv)")
        );

    mvPhiInertial_ =
        fv::convectionScheme<Foam::scalar>::New
        (
            aerosol_.mesh(),
            fields_,
            phiInertial_,
            aerosol_.mesh().divScheme("div(mvConv)")
        );

    // Set the continuous phase diffusivities here, which is aerosol model
    // independent

    if (contDiff_->type() != "none")
    {
        forAll(aerosol_.thermo().contSpecies(), j)
        {
            DCont_[j] = contDiff_->D(j);
        }
    }
    else
    {
        forAll(aerosol_.thermo().contSpecies(), j)
        {
            DCont_[j] *= 0.0;
        }
    }
}

void driftFluxModel::correctPhaseVelocities()
{
    Uc_ = aerosol_.U() - aerosol_.thermo().sumZ()*Ur_;
    Ud_ = aerosol_.U() + aerosol_.thermo().sumY()*Ur_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
