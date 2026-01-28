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

#include "fvCFD.H"
#include "aerosolModel.H"
#include "aerosolThermo.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    autoPtr<aerosolModel> aerosol(aerosolModel::New(mesh));

    aerosolThermo& thermo = aerosol->thermo();

    thermo.correctThermo();
    thermo.correct();

    rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();
    rhoAerosolPhaseThermo& thermoDisp = thermo.thermoDisp();

    basicMultiComponentMixture& compCont = thermoCont.composition();
    basicMultiComponentMixture& compDisp = thermoDisp.composition();

    PtrList<volScalarField>& Y = compCont.Y();
    PtrList<volScalarField>& Z = compDisp.Y();

    Info<< endl << endl;

    Info<< "Continuous species mass fractions:" << nl << endl;

    forAll(mesh.C(), celli)
    {
        Info<< "  Cell " << celli << endl << endl;

        forAll(Y, j)
        {
            Info<< "    " << thermo.species()[j] << " = " << Y[j][celli]
                << endl;
        }
    }

    Info<< endl;

    Info<< "Dispersed species mass fractions:" << nl << endl;

    forAll(mesh.C(), celli)
    {
        Info<< "  Cell " << celli << endl << endl;

        forAll(Z, j)
        {
            Info<< "    " << thermo.species()[j] << " = " << Z[j][celli]
                << endl;
        }
    }

    Info<< endl;

    const scalarField mu(thermo.mu());
    const scalarField rho(thermo.rho());
    const scalarField kappa(thermo.kappa());
    const scalarField alpha(thermo.alpha());

    Info<< "Mixture properties:" << nl << endl;

    forAll(mesh.C(), celli)
    {
        Info<< "  Cell " << celli << endl << endl;

        Info<< "    rho = " << rho[celli] << endl;
        Info<< "    mu = " << mu[celli] << endl;
        Info<< "    kappa = " << kappa[celli] << endl;
        Info<< "    alpha = " << alpha[celli] << endl;
    }

    Info<< endl;

    const scalarField muCont(thermoCont.mu());
    const scalarField rhoCont(thermoCont.rho());
    const scalarField kappaCont(thermoCont.kappa());
    const scalarField alphaCont(thermoCont.alpha());

    Info<< "Continuous mixture properties:" << nl << endl;

    forAll(mesh.C(), celli)
    {
        Info<< "  Cell " << celli << endl << endl;

        Info<< "    rho = " << rhoCont[celli] << endl;
        Info<< "    mu = " << muCont[celli] << endl;
        Info<< "    kappa = " << kappaCont[celli] << endl;
        Info<< "    alpha = " << alphaCont[celli] << endl;
    }

    Info<< endl;

    const scalarField muDisp(thermoDisp.mu());
    const scalarField rhoDisp(thermoDisp.rho());
    const scalarField kappaDisp(thermoDisp.kappa());
    const scalarField alphaDisp(thermoDisp.alpha());

    Info<< "Dispersed mixture properties:" << nl << endl;

    forAll(mesh.C(), celli)
    {
        Info<< "  Cell " << celli << endl << endl;

        Info<< "    rho = " << rhoDisp[celli] << endl;
        Info<< "    mu = " << muDisp[celli] << endl;
        Info<< "    kappa = " << kappaDisp[celli] << endl;
        Info<< "    alpha = " << alphaDisp[celli] << endl;
    }

    Info<< endl;

    const PtrList<scalarField> rhojCont(thermo.rhoCont(thermo.contSpecies()));
    const PtrList<scalarField> mujCont(thermo.muCont(thermo.contSpecies()));
    const PtrList<scalarField> kappajCont(thermo.kappaCont(thermo.contSpecies()));
    const PtrList<scalarField> alphajCont(thermo.alphaCont(thermo.contSpecies()));

    scalarList WCont(thermo.contSpecies().size(), 0.0);

    forAll(WCont, j)
    {
        WCont[j] = thermo.thermoCont().composition().W(j);
    }

    Info<< "Continuous species properties:" << nl << endl;

    forAll(mesh.C(), celli)
    {
        Info<< "  Cell " << celli << endl << endl;

        forAll(Y, j)
        {
            Info<< "    Species " << thermo.species()[j] << endl << endl;

            Info<< "      M = " << WCont[j] << endl;
            Info<< "      rho = " << rhojCont[j][celli] << endl;
            Info<< "      mu = " << mujCont[j][celli] << endl;
            Info<< "      kappa = " << kappajCont[j][celli] << endl;
            Info<< "      alpha = " << alphajCont[j][celli] << endl;

            Info<< endl;
        }
    }

    const PtrList<scalarField> rhojDisp(thermo.rhoDisp(thermo.dispSpecies()));
    const PtrList<scalarField> mujDisp(thermo.muDisp(thermo.dispSpecies()));
    const PtrList<scalarField> kappajDisp(thermo.kappaDisp(thermo.dispSpecies()));
    const PtrList<scalarField> alphajDisp(thermo.alphaDisp(thermo.dispSpecies()));

    scalarList WDisp(thermo.dispSpecies().size(), 0.0);

    forAll(WDisp, j)
    {
        WDisp[j] = thermo.thermoDisp().composition().W(j);
    }

    Info<< "Dispersed species properties:" << nl << endl;

    forAll(mesh.C(), celli)
    {
        Info<< "  Cell " << celli << endl << endl;

        forAll(Z, j)
        {
            Info<< "    Species " << thermo.species()[j] << endl << endl;

            Info<< "      M = " << WDisp[j] << endl;
            Info<< "      rho = " << rhojDisp[j][celli] << endl;
            Info<< "      mu = " << mujDisp[j][celli] << endl;
            Info<< "      kappa = " << kappajDisp[j][celli] << endl;
            Info<< "      alpha = " << alphajDisp[j][celli] << endl;
        }
    }

    Info<< endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}
