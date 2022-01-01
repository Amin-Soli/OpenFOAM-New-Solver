/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "newSchnerrSauer.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(newSchnerrSauer, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeTwoPhaseMixture,
        newSchnerrSauer,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::newSchnerrSauer::newSchnerrSauer
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    n_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("n")),
    dNuc_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("dNuc")),
    Cc_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cc")),
    Cv_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv")),

    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::newSchnerrSauer::rRb
(
    const volScalarField& limitedAlpha2
) const
{
    return pow
    (
        ((4*constant::mathematical::pi*n_)/3)
       *(limitedAlpha2)/(1.0 + alphaNuc() - limitedAlpha2),
        1.0/3.0
    );
}


Foam::dimensionedScalar
Foam::phaseChangeTwoPhaseMixtures::newSchnerrSauer::alphaNuc() const
{
    dimensionedScalar Vnuc = n_*constant::mathematical::pi*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}


Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::newSchnerrSauer::pCoeff
(
    const volScalarField& p
) const
{
    volScalarField limitedAlpha2(min(max(alpha2_, scalar(0)), scalar(1)));
    volScalarField rho
    (
        limitedAlpha2*rho2() + (scalar(1) - limitedAlpha2)*rho1()
    );

    return
        (3*rho1()*rho2())*sqrt(2/(3*rho2()))
       *rRb(limitedAlpha2)/(rho*sqrt(mag(p - pSat()) + 0.01*pSat()));
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::newSchnerrSauer::mDotAlphal() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha2(min(max(alpha2_, scalar(0)), scalar(1)));

    return Pair<tmp<volScalarField>>
    (
        Cc_*limitedAlpha2*pCoeff*max(p - pSat(), p0_),

        Cv_*(1.0 + alphaNuc() - limitedAlpha2)*pCoeff*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::newSchnerrSauer::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha2(min(max(alpha2_, scalar(0)), scalar(1)));
    volScalarField apCoeff(limitedAlpha2*pCoeff);

    return Pair<tmp<volScalarField>>
    (
        Cc_*(1.0 - limitedAlpha2)*pos(p - pSat())*apCoeff,

        (-Cv_)*(1.0 + alphaNuc() - limitedAlpha2)*neg(p - pSat())*apCoeff
    );
}


void Foam::phaseChangeTwoPhaseMixtures::newSchnerrSauer::correct()
{}


bool Foam::phaseChangeTwoPhaseMixtures::newSchnerrSauer::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");

        phaseChangeTwoPhaseMixtureCoeffs_.lookup("n") >> n_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("dNuc") >> dNuc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cc") >> Cc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv") >> Cv_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
