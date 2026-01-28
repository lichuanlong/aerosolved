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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::volVectorField& Foam::driftFluxModel::Ur()
{
    return Ur_;
}

inline const Foam::volVectorField& Foam::driftFluxModel::Ur() const
{
    return Ur_;
}

inline const Foam::volVectorField& Foam::driftFluxModel::Uc() const
{
    return Uc_;
}

inline const Foam::volVectorField& Foam::driftFluxModel::Ud() const
{
    return Ud_;
}

inline const Foam::fv::convectionScheme<Foam::scalar>&
Foam::driftFluxModel::mvPhi() const
{
    return mvPhi_();
}

inline const Foam::fv::convectionScheme<Foam::scalar>&
Foam::driftFluxModel::mvPhiCorr() const
{
    return mvPhiCorr_();
}

inline const Foam::fv::convectionScheme<Foam::scalar>&
Foam::driftFluxModel::mvPhiInertial() const
{
    return mvPhiInertial_();
}

inline
Foam::multivariateSurfaceInterpolationScheme<Foam::scalar>::fieldTable&
Foam::driftFluxModel::fields()
{
    return fields_;
}

inline const
Foam::multivariateSurfaceInterpolationScheme<Foam::scalar>::fieldTable&
Foam::driftFluxModel::fields() const
{
    return fields_;
}

inline Foam::surfaceScalarField& Foam::driftFluxModel::phiCorr()
{
    return phiCorr_;
}

inline const Foam::surfaceScalarField& Foam::driftFluxModel::phiCorr() const
{
    return phiCorr_;
}

inline Foam::surfaceScalarField& Foam::driftFluxModel::phiCorrDiff()
{
    return phiCorrDiff_;
}

inline const Foam::surfaceScalarField& Foam::driftFluxModel::phiCorrDiff() const
{
    return phiCorrDiff_;
}

inline Foam::surfaceScalarField& Foam::driftFluxModel::phiInertial()
{
    return phiInertial_;
}

inline const Foam::surfaceScalarField& Foam::driftFluxModel::phiInertial() const
{
    return phiInertial_;
}

inline Foam::volScalarField& Foam::driftFluxModel::DDisp()
{
    return DDisp_;
}

inline const Foam::volScalarField& Foam::driftFluxModel::DDisp() const
{
    return DDisp_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::driftFluxModel::DCont()
{
    return DCont_;
}

inline const Foam::PtrList<Foam::volScalarField>&
Foam::driftFluxModel::DCont() const
{
    return DCont_;
}

// ************************************************************************* //
