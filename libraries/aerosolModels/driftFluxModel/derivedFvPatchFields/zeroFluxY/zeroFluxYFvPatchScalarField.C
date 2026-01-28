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
    GNU General Public License for more details.s

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/


#include "zeroFluxYFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "fixedSectionalSystem.H"
#include <iostream>
#include <string>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::zeroFluxYFvPatchScalarField::
zeroFluxYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::zeroFluxYFvPatchScalarField::
zeroFluxYFvPatchScalarField
(
    const zeroFluxYFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::zeroFluxYFvPatchScalarField::
zeroFluxYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::zeroFluxYFvPatchScalarField::
zeroFluxYFvPatchScalarField
(
    const zeroFluxYFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::zeroFluxYFvPatchScalarField::
zeroFluxYFvPatchScalarField
(
    const zeroFluxYFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// Code Slice

Foam::label
Foam::zeroFluxYFvPatchScalarField::speciesIndex()
const
{
    // Get AeroSol Model
    if (db().foundObject<aerosolModel>("aerosolProperties"))
    {
        const aerosolModel& aerosol =
              db().lookupObject<aerosolModel>("aerosolProperties");

    	const word suffix(".continuous");
    	const word speciesIndex
    	(
        	this->internalField().name()
    	);

    	forAll(aerosol.thermo().contSpecies(), i)
    	{
        	const word speciesi(aerosol.thermo().contSpecies()[i]);
        	const word speciesIndexi(speciesi+suffix);

		    if (speciesIndex == speciesIndexi)
            {
                return i;
            }
        }

        FatalErrorInFunction
            << "The species to which this BC's field ("
            << this->internalField().name()
            << ") belongs could not be found."
            << "Should be assigned to a continuous species."
            << exit(FatalError);

        return -1;
    }
    else
    {
       FatalErrorInFunction
                << "No Aerosol Model Found"
                << exit(FatalError);

    	return -1;
    }
}

void Foam::zeroFluxYFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get AeroSol Model
    if (db().foundObject<aerosolModel>("aerosolProperties"))
    {
        // Get Aerosol Model
        const aerosolModel& aerosol =
            db().lookupObject<aerosolModel>("aerosolProperties");

        // Mixture Density
        const scalarField& rho =
            aerosol.rho().boundaryField()[patch().index()];

        // Continuous Diffusion Coeffiecient
        const scalarField& D =
            aerosol.drift().DCont()[speciesIndex()].boundaryField()[patch().index()];

        // Mixture Velocity, normal component
        const scalarField Uf
        (
            aerosol.U().boundaryField()[patch().index()] & patch().nf()
        );

        // PhiCorrective
        const scalarField& phiCorr =
            aerosol.drift().phiCorr().boundaryField()[patch().index()];

        // Vcorrective, normal component
        const scalarField Vcorr
        (
            phiCorr/(rho * this->patch().magSf())
        );

        // Molecular Weight
        const volScalarField W
        (
            aerosol.thermo().thermoCont().WMix()
        );

        const surfaceScalarField  W_Diff
        (
            fvc::interpolate(scalar(1.0)/W)*fvc::snGrad(W)
        );

        const scalarField    dW
        (
            W_Diff.boundaryField()[patch().index()]
        );

        // cell size
        const scalarField delta = -patch().delta() & patch().nf();

        // Update s.t. we have 0 deposition at the patch
        operator==
            (patchInternalField()/(scalar(1.0) - delta*dW + delta*(Uf-Vcorr)/(D + scalar(SMALL))));
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}





void Foam::zeroFluxYFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        zeroFluxYFvPatchScalarField
    );
}

// ************************************************************************* //
