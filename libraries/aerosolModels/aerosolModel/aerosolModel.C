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

#include "aerosolModel.H"
#include "aerosolModelGitInfo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aerosolModel, 0);
    defineRunTimeSelectionTable(aerosolModel, dictionary);
}

const Foam::word Foam::aerosolModel::aerosolPropertiesName
(
    "aerosolProperties"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::aerosolModel::versionInfo()
{
    Info<< endl;

    Info<< "###############################################################################" << nl
	    << "####                       Welcome to AeroSolved 2.1                       ####" << nl
	    << "###############################################################################" << nl
	    <<  nl
	    << "Cite: " << nl
	    << "Lucci F., Frederix E.M.A., Kuczaj A.K." << nl
	    << "AeroSolved: Computational fluid dynamics modeling of" << nl
	    << "multispecies aerosol flows with sectional and moment methods," << nl
	    << "Journal of Aerosol Science, Volume 159, (2022)" << nl
	    << "doi:10.1016/j.jaerosci.2021.105854" << nl
	    << nl
	    << "###############################################################################" << nl
	    << nl;

    gitInfo();

	Info<< nl
	    << "###############################################################################"
        << nl << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModel::aerosolModel
(
    const word& modelType,
    const fvMesh& mesh,
    const word& aerosolProperties
)
:
    IOdictionary
    (
        IOobject
        (
            aerosolProperties,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    thermo_(mesh),
    turbulencePtr_(),
    mesh_(mesh),
    coeffs_(modelType == "none" ? *this : subDict(modelType + "Coeffs")),
    modelType_(modelType),
    outputPropertiesPtr_(),
    condensation_(),
    nucleation_(),
    coalescence_(),
    drift_(),
    dMin_
    (
        modelType == "none"
      ? 0.0
      : readScalar(subDict("diameter").lookup("min"))
    ),
    dMax_
    (
        modelType == "none"
      ? 0.0
      : readScalar(subDict("diameter").lookup("max"))
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffs_.lookupOrDefault<scalar>("residualAlpha", 1E-12)
    )
{
    versionInfo();

    read();

    if (!outputPropertiesPtr_.valid())
    {
        const fileName uniformPath(word("uniform")/"aerosolModels");

        outputPropertiesPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "outputProperties",
                    mesh_.time().timeName(),
                    uniformPath,
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }

    if (modelType != "none")
    {
        condensation_ =
            condensationModel::New
            (
                *this,
                subDict("submodels").subDict("condensation")
            );

        nucleation_ =
            nucleationModel::New
            (
                *this,
                subDict("submodels").subDict("nucleation")
            );

        coalescence_ =
            coalescenceModel::New
            (
                *this,
                subDict("submodels").subDict("coalescence")
            );
    }

    drift_.reset(new driftFluxModel(*this, subDict("submodels")));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModel::~aerosolModel()
{
    if (turbulencePtr_)
    {
        turbulencePtr_ = 0;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::aerosolModel::read()
{
    if (regIOobject::read())
    {
        coeffs_ = optionalSubDict(modelType_ + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}

Foam::tmp<Foam::scalarField> Foam::aerosolModel::getRDeltaT()
{
    if(mesh_.foundObject<volScalarField>("rDeltaT"))
    {
        return tmp<scalarField>
        (
            mesh_.lookupObject<volScalarField>("rDeltaT")
        );
    }
    else
    {
        const scalarField& rho = this->rho().field();

        if(!rDeltaT_.valid())
        {
            rDeltaT_.reset(new scalarField(rho.size()));
        }

        rDeltaT_() = 1.0/mesh_.time().deltaTValue();

        return tmp<scalarField>(rDeltaT_());
    }
}

// ************************************************************************* //
