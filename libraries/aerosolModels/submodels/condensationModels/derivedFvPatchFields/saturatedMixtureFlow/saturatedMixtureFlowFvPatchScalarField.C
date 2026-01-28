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

#include "saturatedMixtureFlowFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::saturatedMixtureFlowFvPatchScalarField::solve
(
    PtrList<scalarField>& X,
    PtrList<scalarField>& W,
    const PtrList<scalarField>& S,
    const label facei
) const
{
    // Lists of active species

    scalarList Xa(W.size());
    scalarList Wa(W.size());
    scalarList Sa(W.size());

    forAll(Xa, j)
    {
        Xa[j] = X[j][facei];
        Wa[j] = W[j][facei];
        Sa[j] = S[j][facei];
    }

    // List of inactive species

    scalarList Xi(X.size()-W.size());

    for (label j = W.size(); j < X.size(); j++)
    {
        Xi[j-W.size()] = X[j][facei];
    }

    scalarList XaWa(Xa+Wa);

    bool converged = false;

    for (label iter = 0; iter < maxIter_; iter++)
    {
        const scalar SumX(sum(Xa)+sum(Xi));
        const scalar SumW(sum(Wa));

        const scalarList dXa
        (
            XaWa*Sa*SumX/(SumW+SumX*Sa) - Xa
        );

        if (sum(mag(dXa))/scalar(Xa.size()) < tolerance_)
        {
            converged = true;
            break;
        }
        else
        {
            Xa = Xa + dXa;
            Wa = Wa - dXa;
        }
    }

    // Update solution

    forAll(Xa, j)
    {
        X[j][facei] = Xa[j];
        W[j][facei] = Wa[j];
    }

    return converged;
}

void Foam::saturatedMixtureFlowFvPatchScalarField::checkFields()
{
    if
    (
        !checkFields_
     && db().foundObject<aerosolThermo>("thermophysicalProperties")
    )
    {
        aerosolThermo& thermo =
            db().lookupObjectRef<aerosolThermo>("thermophysicalProperties");

        // Check if all species have this BC

        const PtrList<volScalarField>& Y = thermo.Y();
        const PtrList<volScalarField>& Z = thermo.Z();

        forAll(Y, j)
        {
            if (Y[j].boundaryField()[patch().index()].type() != this->type())
            {
                FatalErrorInFunction
                    << "Boundary condition of " << Y[j].name()
                    << " is not of type " << this->type()
                    << endl << abort(FatalError);
            }
        }

        forAll(Z, j)
        {
            if (Z[j].boundaryField()[patch().index()].type() != this->type())
            {
                FatalErrorInFunction
                    << "Boundary condition of " << Z[j].name()
                    << " is not of type " << this->type()
                    << endl << abort(FatalError);
            }
        }

        // Check if volume fractions are set appropriately

        const speciesTable& species = thermo.species();
        const speciesTable& activeSpecies = thermo.activeSpecies();

        scalarList alphaCont(species.size(), 0);

        forAll(species, i)
        {
            const word specie = species[i];

            alphaCont[i] = alphaCont_.found(specie) ? alphaCont_[specie] : 0;
        }

        if (sum(alphaCont) != 1.0)
        {
            FatalErrorInFunction
                << "Specified continuous phase volume fractions do not add up "
                << "to unity." << nl
                << alphaCont_ << endl << abort(FatalError);
        }

        if (quality_ > 0)
        {
            scalarList alphaDisp(activeSpecies.size(), 0);

            forAll(activeSpecies, i)
            {
                const word specie = activeSpecies[i];

                alphaDisp[i] = alphaDisp_.found(specie) ? alphaDisp_[specie] : 0;
            }

            if (sum(alphaDisp) != 1.0)
            {
                FatalErrorInFunction
                    << "Specified disersed phase volume fractions do not add up "
                    << "to unity." << nl
                    << alphaDisp_ << endl << abort(FatalError);
            }
        }

        checkFields_ = true;
    }


}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturatedMixtureFlowFvPatchScalarField::
saturatedMixtureFlowFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    tolerance_(1e-8),
    maxIter_(1000),
    quality_(0),
    alphaCont_(),
    alphaDisp_(),
    checkFields_(false)
{}


Foam::saturatedMixtureFlowFvPatchScalarField::
saturatedMixtureFlowFvPatchScalarField
(
    const saturatedMixtureFlowFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    tolerance_(ptf.tolerance_),
    maxIter_(ptf.maxIter_),
    quality_(ptf.quality_),
    alphaCont_(ptf.alphaCont_),
    alphaDisp_(ptf.alphaDisp_),
    gamma_(ptf.gamma_),
    checkFields_(false)
{}


Foam::saturatedMixtureFlowFvPatchScalarField::
saturatedMixtureFlowFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-8)),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 1000)),
    quality_(readScalar(dict.lookup("quality"))),
    alphaCont_(),
    alphaDisp_(),
    gamma_(),
    checkFields_(false)
{
    if (quality_ < 0.0 || quality_ >= 1.0)
    {
        FatalErrorInFunction
            << "Flow quality cannot be negative of exceed unity."
            << endl << abort(FatalError);
    }

    // Add continuous phase-specific volume fractions

    if (dict.found("alphaCont"))
    {
        const List<keyType> keys(dict.subDict("alphaCont").keys());

        if (keys.size() < 1)
        {
            FatalErrorInFunction
                << "Continuous phase volume fractions should contain "
                << "at least one species." << endl
                << abort(FatalError);
        }

        forAll(keys, i)
        {
            alphaCont_.insert
            (
                keys[i],
                readScalar(dict.subDict("alphaCont").lookup(keys[i]))
            );
        }
    }
    else
    {
        FatalErrorInFunction
            << "Could not find the volume fractions of the continuous "
            << "phase. " << endl << abort(FatalError);
    }

    // At zero quality, volume fraction does not need to be defined

    if (dict.found("alphaDisp"))
    {
        const List<keyType> keys(dict.subDict("alphaDisp").keys());

        forAll(keys, i)
        {
            alphaDisp_.insert
            (
                keys[i],
                readScalar(dict.subDict("alphaDisp").lookup(keys[i]))
            );
        }
    }
    else if (quality_ > 0)
    {
        FatalErrorInFunction
            << "Could not find the volume fractions of the dispersed "
            << "phase. For non-zero quality at least one specie "
            << "should be defined." << endl << abort(FatalError);
    }

    if (dict.found("gamma"))
    {
        const List<keyType> keys(dict.subDict("gamma").keys());

        forAll(keys, i)
        {
            gamma_.insert
            (
                keys[i],
                readScalar(dict.subDict("gamma").lookup(keys[i]))
            );
        }
    }

    if (!dict.found("value"))
    {
        updateCoeffs();
    }
    else
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
}


Foam::saturatedMixtureFlowFvPatchScalarField::
saturatedMixtureFlowFvPatchScalarField
(
    const saturatedMixtureFlowFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    tolerance_(tppsf.tolerance_),
    maxIter_(tppsf.maxIter_),
    quality_(tppsf.quality_),
    alphaCont_(tppsf.alphaCont_),
    alphaDisp_(tppsf.alphaDisp_),
    gamma_(tppsf.gamma_),
    checkFields_(false)
{}


Foam::saturatedMixtureFlowFvPatchScalarField::
saturatedMixtureFlowFvPatchScalarField
(
    const saturatedMixtureFlowFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    tolerance_(tppsf.tolerance_),
    maxIter_(tppsf.maxIter_),
    quality_(tppsf.quality_),
    alphaCont_(tppsf.alphaCont_),
    alphaDisp_(tppsf.alphaDisp_),
    gamma_(tppsf.gamma_),
    checkFields_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::saturatedMixtureFlowFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    checkFields();

    if (db().foundObject<aerosolThermo>("thermophysicalProperties"))
    {
        aerosolThermo& thermo =
            db().lookupObjectRef<aerosolThermo>("thermophysicalProperties");

        rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();
        rhoAerosolPhaseThermo& thermoDisp = thermo.thermoDisp();

        const basicSpecieMixture& compCont = thermoCont.composition();
        const basicSpecieMixture& compDisp = thermoDisp.composition();

        const speciesTable& species = thermo.species();
        const speciesTable& activeSpecies = thermo.activeSpecies();

        scalarList alphaCont(species.size(), 0);
        scalarList alphaDisp(activeSpecies.size(), 0);

        forAll(species, i)
        {
            const word specie = species[i];

            alphaCont[i] = alphaCont_.found(specie) ? alphaCont_[specie] : 0;
        }

        if (quality_ > 0)
        {
            forAll(activeSpecies, i)
            {
                const word specie = activeSpecies[i];

                alphaDisp[i] = alphaDisp_.found(specie) ? alphaDisp_[specie] : 0;
            }
        }

        const scalarField& p = thermo.p().boundaryField()[patch().index()];
        const scalarField& T = thermo.T().boundaryField()[patch().index()];

        // Mole weights

        scalarList M(species.size(), 0.0);

        forAll(species, j)
        {
            M[j] = thermoCont.composition().W(j);
        }

        // Set activity coefficients. Default value is unity.

        scalarList gamma(activeSpecies.size(), 1.0);

        forAll(activeSpecies, j)
        {
            const word specie = activeSpecies[j];

            if (gamma_.found(specie))
            {
                gamma[j] = gamma_[specie];
            }
        }

        // Saturations

        PtrList<scalarField> S(activeSpecies.size());

        forAll(activeSpecies, j)
        {
            S.set
            (
                j,
                new scalarField
                (
                    gamma[j]
                  * thermoCont.property(activeSpecies[j], "pSat").value(T)
                  / p
                )
            );
        }

        // Mass densities

        PtrList<scalarField> rhoCont(species.size());

        forAll(species, j)
        {
            rhoCont.set(j, new scalarField(patch().size()));

            scalarField& rho = rhoCont[j];

            forAll(rho, facei)
            {
                rho[facei] = compCont.rho(j, p[facei], T[facei]);
            }
        }

        PtrList<scalarField> rhoDisp(activeSpecies.size());

        forAll(activeSpecies, j)
        {
            rhoDisp.set(j, new scalarField(patch().size()));

            scalarField& rho = rhoDisp[j];

            forAll(rho, facei)
            {
                rho[facei] = compDisp.rho(j, p[facei], T[facei]);
            }
        }

        // Initial continuous phase-specific mass fractions

        scalarField sumAlphaRhoCont(patch().size(),0);

        PtrList<scalarField> y(species.size());

        forAll(species, j)
        {
            y.set(j, new scalarField(alphaCont[j]*rhoCont[j]));

            sumAlphaRhoCont += y[j];
        }

        forAll(species, j)
        {
            y[j] /= sumAlphaRhoCont;
        }

        // Initial dispersed phase-specific mass fractions

        scalarField sumAlphaRhoDisp(patch().size(),0);

        PtrList<scalarField> z(activeSpecies.size());

        forAll(activeSpecies, j)
        {
            z.set(j, new scalarField(alphaDisp[j]*rhoDisp[j]));

            sumAlphaRhoDisp += z[j];
        }

        forAll(activeSpecies, j)
        {
            z[j] /= max(sumAlphaRhoDisp,SMALL);
        }

        // Quality in moles

        scalarField moleQuality(patch().size(), 0);

        if (quality_ > 0)
        {
            scalarField MCont(patch().size(),0);
            scalarField MDisp(patch().size(),0);

            forAll(y, j)
            {
                MCont += y[j]/M[j];
            }

            forAll(z, j)
            {
                MDisp += z[j]/M[j];
            }


            MCont = 1.0/MCont;
            MDisp = 1.0/MDisp;

            moleQuality =
                quality_/MDisp
             / (quality_/MDisp + (1.0-quality_)/MCont);
        }

        // Initial continuous mole fractions

        scalarField sumyOverM(patch().size(),0);
        PtrList<scalarField> X(species.size());

        forAll(species, j)
        {
            X.set(j, new scalarField(y[j]/M[j]));

            sumyOverM += y[j]/M[j];
        }

        forAll(species, j)
        {
            X[j] /= sumyOverM;
            X[j] *= (1.0-moleQuality);
        }

        // Initial dispersed mole fractions

        scalarField sumzOverM(patch().size(),0);
        PtrList<scalarField> W(activeSpecies.size());

        forAll(activeSpecies, j)
        {
            W.set(j, new scalarField(z[j]/M[j]));

            sumzOverM += z[j]/M[j];
        }

        forAll(activeSpecies, j)
        {
            W[j] /= max(sumzOverM,SMALL);
            W[j] *= moleQuality;
        }

        // Solve for updated mole fractions

        bool warn = false;

        forAll(*this, facei)
        {
            warn = warn || !this->solve(X,W,S,facei);
        }

        if (warn)
        {
            WarningInFunction
                << "Solution for some faces did not converge"
                << endl;
        }

        // Update the mass fraction of this patch

        const word thisSpecieName = internalField().member();
        const label thisSpecieNum = species[thisSpecieName];

        scalarField sumXWM(patch().size(), 0);

        forAll(X, j)
        {
            sumXWM += X[j]*M[j];
        }

        forAll(W, j)
        {
            sumXWM += W[j]*M[j];
        }

        if (internalField().group() == thermo.phaseNameCont())
        {
            operator==(X[thisSpecieNum]*M[thisSpecieNum]/sumXWM);
        }
        else
        {
            operator==(W[thisSpecieNum]*M[thisSpecieNum]/sumXWM);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::saturatedMixtureFlowFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntry("value", os);

    os.beginBlock("alphaCont");

    forAll(alphaCont_.toc(), i)
    {
        const word key(alphaCont_.toc()[i]);

        os.writeKeyword(key)
            << token::SPACE << alphaCont_[key] << token::END_STATEMENT << nl;
    }

    os.endBlock();

    os.beginBlock("alphaDisp");

    forAll(alphaDisp_.toc(), i)
    {
        const word key(alphaDisp_.toc()[i]);

        os.writeKeyword(key)
            << token::SPACE << alphaDisp_[key] << token::END_STATEMENT << nl;
    }

    os.endBlock();

    os.beginBlock("gamma");

    forAll(gamma_.toc(), i)
    {
        const word key(gamma_.toc()[i]);

        os.writeKeyword(key)
            << token::SPACE << gamma_[key] << token::END_STATEMENT << nl;
    }

    os.endBlock();

    os.writeEntryIfDifferent<scalar>("tolerance", 1e-8, tolerance_);
    os.writeEntryIfDifferent<label>("maxIter", 1000, maxIter_);

    os.writeEntry<scalar>("quality", quality_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        saturatedMixtureFlowFvPatchScalarField
    );
}

// ************************************************************************* //
