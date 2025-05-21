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

#include "saturatedMixtureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::saturatedMixtureFvPatchScalarField::checkFields()
{
    if (!checkFields_)
    {
        const aerosolThermo& thermo =
            db().lookupObject<aerosolThermo>("thermophysicalProperties");

        const speciesTable& contSpecies = thermo.contSpecies();
        const speciesTable& activeSpecies = thermo.activeSpecies();

        specieType_.setSize(contSpecies.size());

        N0_ = 0;
        N1_ = 0;
        N2_ = 0;

        forAll(contSpecies, j)
        {
            const word specie = contSpecies[j];

            if
            (
                thermo.Y()[j].boundaryField()[patch().index()].type()
             == type()
            )
            {
                if (activeSpecies.found(specie))
                {
                    specieType_[j] = 0;
                    N0_++;
                }
                else
                {
                    specieType_[j] = 1;
                    N1_++;
                }
            }
            else
            {
                specieType_[j] = 2;
                N2_++;
            }
        }

        if (N0_ == 1)
        {
            // Only one type 0 specie. Reset the f_ table to have just one entry

            word activeSpecie;

            forAll(activeSpecies, j)
            {
                if (specieType_[j] == 0)
                {
                    activeSpecie = activeSpecies[j];
                    break;
                }
            }

            f_.clear();
            f_.insert(activeSpecie, 1.0);
        }
        else
        {
            scalar sumf = 0.0;

            // Check if the type 0 mole fractions are specified

            forAll(activeSpecies, j)
            {
                if (specieType_[j] == 0)
                {
                    const word specie = activeSpecies[j];

                    if (!f_.found(specie))
                    {
                        FatalErrorInFunction
                            << "Could not find f for " << specie
                            << " in the saturatedMixture BC for "
                            << internalField().name() << endl
                            << abort(FatalError);
                    }
                    else
                    {
                        sumf += f_[specie];
                    }
                }
            }

            // Scale the sum to unity

            forAll(activeSpecies, j)
            {
                if (specieType_[j] == 0)
                {
                    f_[activeSpecies[j]] /= sumf;
                }
            }
        }

        if (N1_ == 0)
        {
            FatalErrorInFunction
                << "No inactive species have the saturatedMixture BC. "
                << "There should be at least one, to act as inert specie."
                << endl << abort(FatalError);
        }
        else if (N1_ == 1)
        {
            // Only one type 1 specie. Reset the g_ table to have just one entry

            word inertSpecie;

            forAll(contSpecies, j)
            {
                if (specieType_[j] == 1)
                {
                    inertSpecie = contSpecies[j];
                    break;
                }
            }

            g_.clear();
            g_.insert(inertSpecie, 1.0);
        }
        else
        {
            scalar sumg = 0.0;

            // Check if the type 1 mole fractions are specified

            forAll(contSpecies, j)
            {
                if (specieType_[j] == 1)
                {
                    const word specie = contSpecies[j];

                    if (!g_.found(specie))
                    {
                        FatalErrorInFunction
                            << "Could not find g for " << specie
                            << " in the saturatedMixture BC for "
                            << internalField().name() << endl
                            << abort(FatalError);
                    }
                    else
                    {
                        sumg += g_[specie];
                    }
                }
            }

            // Scale the sum to unity

            forAll(contSpecies, j)
            {
                if (specieType_[j] == 1)
                {
                    g_[contSpecies[j]] /= sumg;
                }
            }
        }

        checkFields_ = true;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    S_(1.0),
    tolerance_(1e-8),
    checkFields_(false),
    T_(-1.0)
{}


Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const saturatedMixtureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    S_(ptf.S_),
    tolerance_(ptf.tolerance_),
    checkFields_(false),
    T_(ptf.T_),
    f_(ptf.f_),
    g_(ptf.g_)
{}


Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    S_(dict.lookupOrDefault<scalar>("S", 1.0)),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-8)),
    checkFields_(false),
    T_(dict.lookupOrDefault<scalar>("T", -1))
{
    if (dict.found("f"))
    {
        const List<keyType> fKeys(dict.subDict("f").keys());

        forAll(fKeys, i)
        {
            f_.insert
            (
                fKeys[i],
                readScalar(dict.subDict("f").lookup(fKeys[i]))
            );
        }
    }

    if (dict.found("g"))
    {
        const List<keyType> gKeys(dict.subDict("g").keys());

        forAll(gKeys, i)
        {
            g_.insert
            (
                gKeys[i],
                readScalar(dict.subDict("g").lookup(gKeys[i]))
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


Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const saturatedMixtureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    S_(tppsf.S_),
    tolerance_(tppsf.tolerance_),
    checkFields_(false),
    T_(tppsf.T_),
    f_(tppsf.f_),
    g_(tppsf.g_)
{}


Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const saturatedMixtureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    S_(tppsf.S_),
    tolerance_(tppsf.tolerance_),
    checkFields_(false),
    T_(tppsf.T_),
    f_(tppsf.f_),
    g_(tppsf.g_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::saturatedMixtureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    checkFields();

    if
    (
        patch().size() > 0
     && db().foundObject<aerosolThermo>("thermophysicalProperties")
    )
    {
        aerosolThermo& thermo =
            db().lookupObjectRef<aerosolThermo>("thermophysicalProperties");

        rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();

        const speciesTable& contSpecies = thermo.contSpecies();
        const speciesTable& dispSpecies = thermo.dispSpecies();

        const scalarField& p = thermo.p().boundaryField()[patch().index()];
        const scalarField& Tp = thermo.T().boundaryField()[patch().index()];

        const scalarField T
        (
            T_ > 0.0
          ? scalarField(patch().size(), T_)
          : Tp
        );

        // Set mole weights

        scalarList W(contSpecies.size(), 0.0);

        forAll(contSpecies, j)
        {
            W[j] = thermoCont.composition().W(j);
        }

        // Set type 0 and type 2 specie saturation pressures. Leave unset for
        // type 1 species

        PtrList<scalarField> pSat(contSpecies.size());

        forAll(contSpecies, j)
        {
            const word specie = contSpecies[j];

            if (specieType_[j] == 0 || specieType_[j] == 2)
            {
                pSat.set
                (
                    j,
                    new scalarField
                    (
                        thermoCont.property(specie, "pSat").value(T)
                    )
                );
            }
        }

        // Set the activity coefficient fields

        PtrList<scalarField> w(dispSpecies.size());

        forAll(dispSpecies, j)
        {
            w.set(j, new scalarField(patch().size(), f_[dispSpecies[j]]));
        }

        PtrList<scalarField> gamma
        (
            thermo.activity().gamma
            (
                dispSpecies,
                p,
                Tp,
                w
            )
        );

        // Compute the sum of continuous mass fractions as one minus the sum of
        // dispersed mass fractions, from the interior cells. This implicitly
        // applies a zero-gradient BC to the dispersed mass fraction fields.

        scalarField sumY(patch().size(), 1.0);

        forAll(dispSpecies, j)
        {
            sumY -=
                thermo.Z()[j]
               .boundaryField()[patch().index()]
               .patchInternalField();
        }

        // Type 2 mass fractions w.r.t. the continuous phase

        PtrList<scalarField> y2(N2_);

        label k = 0;

        forAll(contSpecies, j)
        {
            if (specieType_[j] == 2)
            {
                y2.set
                (
                    k++,
                    new scalarField
                    (
                        thermo.Y()[j].boundaryField()[patch().index()]
                      / sumY
                    )
                );
            }
        }

        // Initialize a list of continuous phase mole fraction fields for all
        // species

        PtrList<scalarField> x(contSpecies.size());

        forAll(contSpecies, j)
        {
            x.set(j, new scalarField(patch().size(), 0.0));
        }

        // Initial guess of type 2 mole fractions w.r.t. the continuous phase

        k = 0;

        forAll(contSpecies, j)
        {
            if (specieType_[j] == 2)
            {
                x[j] = y2[k++];
            }
        }

        // Iterate on the computation of the continuous phase composition and
        // associated mean molar mass M, which is used in the convergence
        // criterion

        scalarField M(patch().size(), 0.0);

        label NIter = 100;

        for (label iter = 0; iter < NIter; iter++)
        {
            // Compute liquid pool type 2 mole fractions from Raoult's law and
            // S0

            PtrList<scalarField> w2(N2_);
            scalarField S0(patch().size(), S_);

            k = 0;

            forAll(contSpecies, j)
            {
                if (specieType_[j] == 2)
                {
                    w2.set(k, new scalarField(x[j]*gamma[j]*p/pSat[j]));
                    S0 -= w2[k];
                    k++;
                }
            }

            // Limit S0 so that it is great than or equal to zero, to avoid
            // negative mole fractions in the pool

            S0 = max(S0,0.0);

            // Compute liquid pool type 0 mole fractions from f

            PtrList<scalarField> w0(N0_);

            k = 0;

            forAll(contSpecies, j)
            {
                if (specieType_[j] == 0)
                {
                    w0.set(k++, new scalarField(S0*f_[contSpecies[j]]));
                }
            }

            // From the liquid pool mole fractions we can now compute the
            // continuous phase mass fractions using Raoult's law again

            k = 0;

            forAll(contSpecies, j)
            {
                if (specieType_[j] == 0)
                {
                    x[j] = w0[k++]*gamma[j]*pSat[j]/p;
                }
            }

            // Compute the total mole fraction of type 1 species in the
            // continuous phase

            scalarField T1(patch().size(), 1.0);

            forAll(contSpecies, j)
            {
                if (specieType_[j] == 0 || specieType_[j] == 2)
                {
                    T1 -= x[j];
                }
            }

            // Limit T1 to avoid negative mole fractions. Note that if this
            // limit is reached, the mixture is probably not very realistic
            // because we will have no type 1 (inert) species in it

            T1 = max(T1,0.0);

            // Compute the type 1 mole fractions in the continuous phase, from g

            forAll(contSpecies, j)
            {
                if (specieType_[j] == 1)
                {
                    x[j] = T1*g_[contSpecies[j]];
                }
            }

            // Update the continuous phase mixture mean molar mass

            scalarField MNew(patch().size(), 0.0);

            forAll(contSpecies, j)
            {
                MNew += W[j]*x[j];
            }


            if (average(mag(MNew - M)) < tolerance_ && iter > 2)
            {
                break;
            }
            else if (iter == NIter-1)
            {
                WarningInFunction
                    << "Calculation of saturated mixture composition failed "
                    << "to converge. You may try to change the tolerance, "
                    << "which is currently set to " << tolerance_ << ". "
                    << "Achieved residual = " << average(mag(MNew - M))
                    << endl;
            }
            else
            {
                M = MNew;

                // Update the type 2 continuous phase mole fractions using the new
                // continuous mixture mean molar mass, for the next iteration

                k = 0;

                forAll(contSpecies, j)
                {
                    if (specieType_[j] == 2)
                    {
                        x[j] = y2[k++]/W[j]*M;
                    }
                }
            }

        }
        // Set the mass fraction w.r.t. the total mixture of this specie

        const word thisSpecieName = internalField().member();
        const label thisSpecieNum = contSpecies[thisSpecieName];

        operator==(x[thisSpecieNum]*W[thisSpecieNum]/M*sumY);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::saturatedMixtureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntry("value", os);

    os.beginBlock("f");

    forAll(f_.toc(), i)
    {
        const word key(f_.toc()[i]);

        os.writeKeyword(key)
            << token::SPACE << f_[key] << token::END_STATEMENT << nl;
    }

    os.endBlock();

    os.beginBlock("g");

    forAll(g_.toc(), i)
    {
        const word key(g_.toc()[i]);

        os.writeKeyword(key)
            << token::SPACE << g_[key] << token::END_STATEMENT << nl;
    }

    os.endBlock();

    os.writeEntryIfDifferent<scalar>("S", 1.0, S_);
    os.writeEntryIfDifferent<scalar>("tolerance", 1e-8, tolerance_);
    os.writeEntryIfDifferent<scalar>("T", -1.0, T_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        saturatedMixtureFvPatchScalarField
    );
}

// ************************************************************************* //
