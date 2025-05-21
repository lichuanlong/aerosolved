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

\file setLogNormal.C

\brief Utility to set a saturated mixture of species.

This utility acts much like the saturatedMixture BC, with the difference that
there are no 'type 2' species. All active species are of type 0 and all inactive
species are of type 1. Parameters can be set in the setSaturatedMixtureDict
dictionary in system. For further documentation, see the saturatedMixture BC.

*/

#include "argList.H"
#include "fvCFD.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    aerosolThermo thermo(mesh);

    rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();

    const speciesTable& contSpecies = thermo.contSpecies();
    const speciesTable& activeSpecies = thermo.activeSpecies();
    const speciesTable& inactiveSpecies = thermo.inactiveSpecies();

    IOdictionary dict
    (
        IOobject
        (
            "setSaturatedMixtureDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (activeSpecies.size() > 1 && !dict.headerOk())
    {
        FatalErrorInFunction
            << "There are more than one inactive species, so their mole "
            << "fractions w.r.t. the 'liquid pool mixture' (f) should be "
            << "specified in system/setSaturatedMixtureDict, which cannot "
            << "be found." << endl << abort(FatalError);
    }

    if (inactiveSpecies.size() > 1 && !dict.headerOk())
    {
        FatalErrorInFunction
            << "There are more than one active species, so their mole "
            << "fractions w.r.t. the inert mixture (f) should be specified "
            << "in system/setSaturatedMixtureDict, which cannot be found."
            << endl << abort(FatalError);
    }

    // Read S

    scalar S = 1.0;

    if (dict.headerOk() && dict.found("S"))
    {
        S = readScalar(dict.lookup("S"));
    }

    // Read f if there are more than one active species

    scalarList f(activeSpecies.size(), 1.0);

    if (activeSpecies.size() > 1)
    {
        scalar sumf = 0.0;

        forAll(activeSpecies, j)
        {
            if (dict.subDict("f").found(activeSpecies[j]))
            {
                f[j] = readScalar(dict.subDict("f").lookup(activeSpecies[j]));

                sumf += f[j];
            }
        }

        f = f/sumf;
    }

    // Read g if there are more than one inactive species

    scalarList g(inactiveSpecies.size(), 1.0);

    if (inactiveSpecies.size() > 1)
    {
        scalar sumg = 0.0;

        forAll(inactiveSpecies, j)
        {
            if (dict.subDict("g").found(inactiveSpecies[j]))
            {
                g[j] = readScalar(dict.subDict("g").lookup(inactiveSpecies[j]));

                sumg += g[j];
            }
        }

        g = g/sumg;
    }

    // Read activity coefficients

    scalarList gamma(activeSpecies.size(), 1.0);

    if (dict.found("gamma"))
    {
        forAll(activeSpecies, j)
        {
            if (dict.subDict("gamma").found(activeSpecies[j]))
            {
                gamma[j] =
                    readScalar(dict.subDict("gamma").lookup(activeSpecies[j]));
            }
        }
    }

    // Pressure and temperature field references

    const scalarField& p = thermo.p().field();
    const scalarField& T = thermo.T().field();

    // Set mole weights

    scalarList W(contSpecies.size(), 0.0);

    forAll(contSpecies, j)
    {
        W[j] = thermoCont.composition().W(j);
    }

    // List of continuous phase mole fractions

    PtrList<scalarField> xs(contSpecies.size());

    // Set active species mole fractions w.r.t. the continuous phase, using
    // Raoult's law

    scalarField sumx(p.size(), 0.0);

    forAll(activeSpecies, j)
    {
        xs.set
        (
            j,
            S*f[j]*gamma[j]*thermoCont.property(activeSpecies[j], "pSat").value(T)/p
        );

        sumx += xs[j];
    }

    // Set inactive species mole fractions w.r.t. the continuous phase

    forAll(inactiveSpecies, j)
    {
        xs.set(activeSpecies.size()+j, (1.0-sumx)*g[j]);
    }

    // Compute mixture mean molar mass

    scalarField M(p.size(), 0.0);

    forAll(contSpecies, j)
    {
        M += xs[j]*W[j];
    }

    // Store mass fractions w.r.t. the total mixture

    const scalarField sumZ(thermo.sumZ()().field());

    forAll(contSpecies, j)
    {
        thermo.Y()[j].field() = xs[j]*W[j]/M*(1.0-sumZ);
        thermo.Y()[j].write();
    }

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
