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

#include "dispersedInertialDriftModel.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dispersedInertialDriftModel, 0);
defineRunTimeSelectionTable(dispersedInertialDriftModel, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

volVectorField& dispersedInertialDriftModel::getCachedField
(
    const word sizeName
)
{
    const word name = this->velocityFieldName(sizeName);

    if (!cache_.found(name))
    {
        // Create or read field. Print outcome of this operation so that the
        // user will be aware of possible cache misses. Explicitly read from
        // start time, because we may already be advanced by one deltaT once we
        // reach this point.

        const fvMesh& mesh = aerosol_.mesh();

        const word startTimeName
        (
            mesh.time().timeName
            (
                mesh.time().startTime().value()
            )
        );

        IOobject fieldHeader
        (
            name,
            startTimeName,
            mesh,
            IOobject::NO_READ
        );

        if (fieldHeader.typeHeaderOk<volVectorField>(true))
        {
            Info<< "Reading drift velocity field " << name
                << " and storing in cache." << endl;

            cache_.insert
            (
                name,
                autoPtr<volVectorField>
                (
                    new volVectorField
                    (
                        IOobject
                        (
                            name,
                            startTimeName,
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh
                    )
                )
            );
        }
        else
        {
            Info<< "Creating drift velocity field " << name
                << " and storing in cache." << endl;

            cache_.insert
            (
                name,
                autoPtr<volVectorField>
                (
                    new volVectorField
                    (
                        IOobject
                        (
                            name,
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        VPtr_()
                    )
                )
            );
        }
    }

    return *cache_[name];
}

void dispersedInertialDriftModel::readBaseField()
{
    VPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "V",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            aerosol_.mesh()
        )
    );
}

void dispersedInertialDriftModel::readBaseFieldIfPresent()
{
    VPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "V",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            aerosol_.mesh(),
            dimensionedVector("V", dimVelocity, vector::zero)
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dispersedInertialDriftModel::dispersedInertialDriftModel
(
    const word& modelType,
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    aerosolSubModelBase(aerosol, dict, typeName, modelType),
    VPtr_(),
    cache_(0),
    VMax_(dict.lookupOrDefault<scalar>("VMax", 1e99)),
    g_
    (
        IOobject
        (
            "g",
            aerosol.mesh().time().constant(),
            aerosol.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimVelocity/dimTime, vector::zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dispersedInertialDriftModel::~dispersedInertialDriftModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

word dispersedInertialDriftModel::velocityFieldName(const word sizeName) const
{
    return IOobject::groupName("V", sizeName);
}

tmp<volScalarField>
dispersedInertialDriftModel::Re
(
    const volScalarField& d,
    const volVectorField& V
) const
{
    tmp<volScalarField> tRe
    (
        new volScalarField
        (
            IOobject
            (
                "Re",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aerosol_.mesh(),
            dimensionedScalar("Re", dimless, 0.0)
        )
    );

    const volScalarField& rhog = aerosol_.thermo().thermoCont().rho();
    const volScalarField& mug = aerosol_.thermo().thermoCont().mu();

    tRe.ref() = d*rhog*mag(V)/mug;

    return tRe;
}

void dispersedInertialDriftModel::limit(volVectorField& V) const
{
    const dimensionedScalar VMin("VMin", V.dimensions(), SMALL);
    const dimensionedScalar VMax("VMax", V.dimensions(), VMax_);

    V *= min(mag(V),VMax)/max(mag(V),VMin);
    V.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
