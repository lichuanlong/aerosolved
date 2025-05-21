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

#include "fixedSectionalSystem.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedSectionalSystem::fixedSectionalSystem
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    aerosolSubModelBase
    (
        aerosol,
        dict,
        "fixedSectionalSystem",
        "fixedSectionalSystem"
    ),
    regIOobject
    (
        IOobject
        (
            "fixedSectionalSystem",
            aerosol_.mesh().time().constant()/"fixedSectionalSystem",
            aerosol_.mesh()
        )
    ),
    distribution_(),
    interpolation_(),
    M_
    (
        IOobject
        (
            "M",
            aerosol.mesh().time().timeName(),
            aerosol.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aerosol.mesh(),
        dimensionedScalar("M", dimless/dimMass, 0.0)
    )
{
    distribution_ = sectionalDistribution::New
    (
        aerosol,
        dict.subDict("distribution"),
        dimMass
    );

    forAll(distribution_(), i)
    {
        fields_.add(distribution_()[i].M());
    }

    Info<< endl << "Representative size list [kg] =" << endl << "(" << endl;

    forAll(distribution_->x(), i)
    {
        Info<< token::TAB << distribution_->x()[i] << endl;
    }

    Info<< ")" << endl << endl;

    interpolation_ = sectionalInterpolation::New
    (
        aerosol,
        distribution_(),
        dict.subDict("interpolation")
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fixedSectionalSystem::~fixedSectionalSystem()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> fixedSectionalSystem::d(const label i) const
{
    tmp<volScalarField> td
    (
        new volScalarField
        (
            IOobject
            (
                "d",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aerosol_.mesh(),
            dimensionedScalar("d", dimLength, 0.0)
        )
    );

    volScalarField& d = td.ref();

    const volScalarField& rhol = aerosol_.thermo().thermoDisp().rho();

    forAll(d.field(), celli)
    {
        d.field()[celli] = distribution_()[i].d(rhol.field()[celli]);
    }

    forAll(d.boundaryField(), patchi)
    {
        d.boundaryFieldRef()[patchi] =
            distribution_()[i].d(rhol.boundaryField()[patchi]);
    }

    return td;
}

tmp<scalarField> fixedSectionalSystem::d
(
    const label i,
    const label patchi
) const
{
    const scalarField rhol(aerosol_.thermo().thermoDisp().rho(patchi));

    tmp<scalarField> td
    (
        new scalarField(distribution_()[i].d(rhol))
    );

    return td;
}

tmp<volScalarField> fixedSectionalSystem::meanDiameter
(
    const scalar p,
    const scalar q
) const
{
    // We have to make sure that all boundaries of all M fields are properly
    // initialized. This requires a const cast.

    sectionalDistribution& dist =
        const_cast<sectionalDistribution&>(distribution_());

    dist.correctBoundaryConditions();

    const volScalarField d0(d(0));

    const dimensionedScalar zeroM("M", M_.dimensions(), 0.0);
    const volScalarField& M0 = dist[0].M();

    volScalarField dp(Foam::pow(d0,p)*M0);
    volScalarField dq(Foam::pow(d0,q)*M0);

    for (label i = 1; i < distribution_->size(); i++)
    {
        const volScalarField di(d(i));
        const volScalarField& Mi = dist[i].M();

        dp += Foam::pow(di,p)*Mi;
        dq += Foam::pow(di,q)*Mi;
    }

    dp.max(dimensionedScalar("d", dp.dimensions(), 0.0));
    dq.max(dimensionedScalar("d", dq.dimensions(), SMALL));

    return Foam::pow(dp/dq, 1.0/(p-q));
}

tmp<scalarField> fixedSectionalSystem::meanDiameter
(
    const scalar p,
    const scalar q,
    const label patchi
) const
{
    const sectionalDistribution& dist = distribution_();

    const scalarField d0(d(0,patchi));

    const scalarField& M0 = dist[0].M().boundaryField()[patchi];

    scalarField dp(Foam::pow(d0,p)*M0);
    scalarField dq(Foam::pow(d0,q)*M0);

    for (label i = 1; i < distribution_->size(); i++)
    {
        const scalarField di(d(i,patchi));

        const scalarField& Mi = dist[i].M().boundaryField()[patchi];

        dp += Foam::pow(di,p)*Mi;
        dq += Foam::pow(di,q)*Mi;
    }

    return Foam::pow(max(dp,0.0)/max(dq,SMALL), 1.0/(p-q));
}

tmp<volScalarField> fixedSectionalSystem::medianDiameter
(
    const scalar p
) const
{
    tmp<volScalarField> tpMD
    (
        new volScalarField
        (
            IOobject
            (
                word("medianDiameter("+Foam::name(p)+")"),
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aerosol_.mesh(),
            dimensionedScalar("zero", dimLength, 0.0)
        )
    );

    scalarField& pMD = tpMD.ref();

    const scalarField& rhol = aerosol_.thermo().thermoDisp().rho();

    forAll(pMD, celli)
    {
        const scalar pMX(distribution_().median(celli, p));

        pMD[celli] =
            Foam::pow
            (
                pMX/rhol[celli]*6.0/constant::mathematical::pi,
                1.0/3.0
            );
    }

    return tpMD;
}

tmp<volScalarField> fixedSectionalSystem::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                "sumf",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            distribution_()[0].M()*distribution_()[0].xd()
        )
    );

    volScalarField& alpha = talpha.ref();

    for (label i = 1; i < distribution_->size(); i++)
    {
        alpha += distribution_()[i].M()*distribution_()[i].xd();
    }

    return talpha;
}

tmp<volScalarField> fixedSectionalSystem::sumM() const
{
    tmp<volScalarField> tsumM
    (
        new volScalarField
        (
            IOobject
            (
                "sumM",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            distribution_()[0].M()
        )
    );

    volScalarField& sumM = tsumM.ref();

    for (label i = 1; i < distribution_->size(); i++)
    {
        sumM += distribution_()[i].M();
    }

    return tsumM;
}

void fixedSectionalSystem::rescale()
{
    PtrList<section>& sections = distribution_->sections();

    forAll(sections, i)
    {
        sections[i].M().max(0.0);
    }

    volScalarField alphaFromM(this->alpha());

    alphaFromM.max
    (
        dimensionedScalar("alpha", alphaFromM.dimensions(), sqrt(VSMALL))
    );

    const volScalarField alphaFromZ(aerosol_.thermo().sumZ());

    const scalarField delta(alphaFromM.field()-alphaFromZ.field());

    Info<< "fixedSectionalSystem: mass fraction difference, min, max = "
        << gSum(delta*aerosol_.mesh().V().field())
         / gSum(aerosol_.mesh().V().field())
        << ", " << gMin(delta)
        << ", " << gMax(delta)
        << endl;

    // Correct using a mass weighting, such that most correction is in the
    // sections with most mass. This is of the form:
    //
    // M_{i,new} = M_i*(1 + corr*x_i*M_i)
    //
    // in which the corr factor can be inferred from the requirement that
    //
    // sum_i M_{i,new}*x_i = sum_j Z_j
    //

    volScalarField alphaFromMSqr
    (
        IOobject
        (
            "alphaFromMSqr",
            aerosol_.mesh().time().timeName(),
            aerosol_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aerosol_.mesh(),
        dimensionedScalar("alphaFromMSqr", dimless, 0)
    );

    forAll(sections, i)
    {
        alphaFromMSqr += sqr(sections[i].xd()*sections[i].M());
    }

    alphaFromMSqr.max
    (
        dimensionedScalar("alphaSqr", alphaFromMSqr.dimensions(), VSMALL)
    );

    const volScalarField corr
    (
        (alphaFromZ - alphaFromM)/alphaFromMSqr
    );

    forAll(sections, i)
    {
        sections[i].M() *= (1.0 + corr*sections[i].xd()*sections[i].M());
    }

    // Finally, linearly scale to assure positivity and correct BCs

    forAll(sections, i)
    {
        sections[i].M().max(0.0);
    }

    alphaFromM = this->alpha();

    alphaFromM.max
    (
        dimensionedScalar("alpha", alphaFromM.dimensions(), sqrt(VSMALL))
    );

    forAll(sections, i)
    {
        sections[i].M() *= alphaFromZ/alphaFromM;
        sections[i].M().correctBoundaryConditions();
    }
}

void fixedSectionalSystem::collect()
{
    PtrList<section>& sections = distribution_->sections();

    M_ == dimensionedScalar(M_.dimensions(), Zero);

    forAll(sections, i)
    {
        M_ += sections[i].M();
    }
}

void fixedSectionalSystem::generateCoalescencePairs()
{
    PtrList<section>& sections = distribution_->sections();

    const label P(sections.size());

    coalescencePairs_.clear();
    coalescencePairs_.setSize(label(0.5*P*(P+1)));

    label k(-1);

    forAll(sections, i)
    {
        for (label j = i; j < sections.size(); j++)
        {
            const scalar s(sections[i].x()+sections[j].x());

            k++;

            coalescencePairs_.set
            (
                k,
                new coalescencePair(i, j, s, interpolation_->interp(s))
            );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
