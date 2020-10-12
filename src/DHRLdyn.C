/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

#include "DHRLdyn.H"

#include "addToRunTimeSelectionTable.H"

#include "wallDist.H"

//#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DHRLdyn, 0);
addToRunTimeSelectionTable(LESModel, DHRLdyn, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

dimensionedScalar DHRLdyn::cD
(
    const volSymmTensorField& D
) const
{
    const volSymmTensorField MM
    (
        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D))
    );

    dimensionedScalar MMMM = average(magSqr(MM));

    if (MMMM.value() > VSMALL)
    {
        tmp<volSymmTensorField> LL =
            dev(filter_(sqr(U())) - (sqr(filter_(U()))));

        return 0.5*average(LL && MM)/MMMM;
    }
    else
    {
        return 0.0;
    }
}


dimensionedScalar DHRLdyn::cI
(
    const volSymmTensorField& D
) const
{
    const volScalarField mm
    (
        sqr(delta())*(4*sqr(mag(filter_(D))) - filter_(sqr(mag(D))))
    );

    dimensionedScalar mmmm = average(magSqr(mm));

    if (mmmm.value() > VSMALL)
    {
        tmp<volScalarField> KK =
            0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

        return average(KK*mm)/mmmm;
    }
    else
    {
        return 0.0;
    }
}

tmp<volScalarField> DHRLdyn::F1(const volScalarField& CDkOmega) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


tmp<volScalarField> DHRLdyn::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
DHRLdyn::DHRLdyn
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    ///// -- SMAGORINSKY COEFFS -- /////
     ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),
     ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce",
            coeffDict_,
            1.048
        )
    ),
    ////////////////////////////////////////////
    y_(mesh_),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_

    ),

    ks
    (
        IOobject
        (
            "ks",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        k_

    ),

    // Define a lower bound for omega as SMALL = 1e-15 -- ERIC
    omegaMin_("omegaMin", dimless/dimTime, SMALL),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_

    ),

Uavg
(
    IOobject
    (
        "Uavg",
        runTime_.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U_//Uavg set to U initially
),

tauSgsAvg
(
    IOobject
    (
        "tauSgsAvg",
        runTime_.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedSymmTensor("tauSgsAvg",sqr(dimLength/dimTime), symmTensor::zero)
),

tauRANS
(
    IOobject
    (
        "tauRANS",
        runTime_.timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedSymmTensor("tauRANS",sqr(dimLength/dimTime), symmTensor::zero)
),
RTP
(
    IOobject
    (
        "RTP",
        runTime_.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("RTP",sqr(dimLength)/pow(dimTime,3), pTraits<scalar>::zero)
),
SP
(
    IOobject
    (
        "SP",
        runTime_.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("SP",sqr(dimLength)/pow(dimTime,3), pTraits<scalar>::zero)
),

RP
(
    IOobject
    (
        "RP",
        runTime_.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("RP",sqr(dimLength)/pow(dimTime,3), pTraits<scalar>::zero)
),

alpha
(
    IOobject
    (
        "alpha",
        runTime_.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("alpha", dimless, pTraits<scalar>::zero)
),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_

    ),
nuSgsR
    (
        IOobject
        (
            "nuSgsR",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nuSgs_

    ),
nuSgsL
    (
        IOobject
        (
            "nuSgsL",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nuSgs_

    ),

UPrimeAvg11
    (
    	IOobject
        (
            "UPrimeAvg11",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("UPrimeAvg11", sqr(dimLength/dimTime), pTraits<scalar>::zero)

    ),

UPrimeAvg22
    (
    	IOobject
        (
            "UPrimeAvg22",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("UPrimeAvg22", sqr(dimLength/dimTime), pTraits<scalar>::zero)

    ),

UPrimeAvg33
    (
    	IOobject
        (
            "UPrimeAvg33",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("UPrimeAvg33", sqr(dimLength/dimTime), pTraits<scalar>::zero)

    ),

UPrimeAvg12
    (
    	IOobject
        (
            "UPrimeAvg12",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("UPrimeAvg12", sqr(dimLength/dimTime), pTraits<scalar>::zero)

    ),

UPrimeAvg13
    (
    	IOobject
        (
            "UPrimeAvg13",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("UPrimeAvg13", sqr(dimLength/dimTime), pTraits<scalar>::zero)

    ),

UPrimeAvg23
    (
    	IOobject
        (
            "UPrimeAvg23",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("UPrimeAvg23", sqr(dimLength/dimTime), pTraits<scalar>::zero)

    ),

 filterPtr_(LESfilter::New(U_.mesh(), coeffDict())),
    filter_(filterPtr_())

{
Info << "\nDynamic Hybrid RANS-LES (k-w SST & Dynamic Smagorinsky)\n" << endl;
    omegaMin_.readIfPresent(*this);
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

 nuSgsR =
    (
        a1_*k_
      / max
        (
            a1_*omega_,
            F2()*sqrt(2.0)*mag(symm(fvc::grad(Uavg)))
        )
    );
    nuSgsR.correctBoundaryConditions();


    printCoeffs();


}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> DHRLdyn::B() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "B",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            //Explicit definition of the Reynolds stress tensor
            //this is not used in the momentum equation; rather we have
            //div(B) = divDevBeff (further down) which goes into the solver
            //the function divDevBeff places an implicit treatment on calculating div(B)
            ((2.0/3.0)*I)*((1.0 - alpha)*ks + alpha*k_) -
            (1.0 - alpha)*tauSgsAvg + alpha*tauRANS,
            k_.boundaryField().types()

        )
    );
}

tmp<volSymmTensorField> DHRLdyn::devBeff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoBeff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           //This is needed for computing wall shear stress
           //where tau_wall = (nu + nut)*(du/dy)_wall -- Eric
           - nuEff()*((1 - alpha)*dev(twoSymm(fvc::grad(U_)))
           + alpha*dev(twoSymm(fvc::grad(Uavg))))
        )
    );
}

tmp<fvVectorMatrix> DHRLdyn::divDevBeff(volVectorField& U) const
{
    return
    (

      - fvm::laplacian(nu() + (1-alpha)*nuSgsL, U)
      - fvc::div((nu() + (1-alpha)*nuSgsL)*dev(T(fvc::grad(U))))
      - fvc::div(alpha*nuSgsR*dev(T(fvc::grad(Uavg))))

    );
}


bool DHRLdyn::read()
{
    if (LESModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        ce_.readIfPresent(coeffDict());
        ck_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}

void DHRLdyn::correct()
{
    LESModel::correct();
/*
    if (!turbulence_)
    {
        return;
    }
*/
    if (mesh_.changing())
    {
        y_.correct();
    }

const volSymmTensorField D(dev(symm(fvc::grad(U_))));


scalar time = runTime_.value();
scalar dt = runTime_.deltaT().value();
scalar startAvg = dt;

if (time > startAvg )
{
    if ( (time - dt) < startAvg )
    {
        Uavg = U_;

    }
    else
    {
        Uavg = ( Uavg.oldTime()*(time-startAvg) + U_*dt ) / (time-startAvg + dt);
        //Info << "\nDone calculating averaged velocity field.\n" << endl;

    }
}

    volVectorField UPrime = U_ - Uavg;
    //Info << "\nDone calculating velocity fluctuation.\n" << endl;

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    //Calculate stresses
    volScalarField UPrime11 = UPrime.component(0)*UPrime.component(0);
    volScalarField UPrime22 = UPrime.component(1)*UPrime.component(1);
    volScalarField UPrime33 = UPrime.component(2)*UPrime.component(2);
    volScalarField UPrime12 = UPrime.component(0)*UPrime.component(1);
    volScalarField UPrime13 = UPrime.component(0)*UPrime.component(2);
    volScalarField UPrime23 = UPrime.component(1)*UPrime.component(2);

    //Info << "\nDone calculating stress components.\n" << endl;

    // Calculate the average of the velocity fluctuations //
    if (time > startAvg )
{
    if ( (time - dt) < startAvg )
    {
        UPrimeAvg11 = UPrime11;
        UPrimeAvg22 = UPrime22;
        UPrimeAvg33 = UPrime33;
        UPrimeAvg12 = UPrime12;
        UPrimeAvg13 = UPrime13;
        UPrimeAvg23 = UPrime23;
    }
    else
    {
    UPrimeAvg11 = ( UPrimeAvg11.oldTime()*(time-startAvg) + UPrime11*dt ) / ((time-startAvg) + dt);

    UPrimeAvg22 = ( UPrimeAvg22.oldTime()*(time-startAvg) + UPrime22*dt ) / ((time-startAvg) + dt);


    UPrimeAvg33 = ( UPrimeAvg33.oldTime()*(time-startAvg) + UPrime33*dt ) / ((time-startAvg) + dt);

    UPrimeAvg12 = ( UPrimeAvg12.oldTime()*(time-startAvg) + UPrime12*dt ) / ((time-startAvg) + dt);

    UPrimeAvg13 = ( UPrimeAvg13.oldTime()*(time-startAvg) + UPrime13*dt ) / ((time-startAvg) + dt);

    UPrimeAvg23 = ( UPrimeAvg23.oldTime()*(time-startAvg) + UPrime23*dt ) / ((time-startAvg) + dt);

    //Info << "\nDone composing UPrime.\n" << endl;
    }
}

volSymmTensorField UPrimeAvg = volSymmTensor(UPrimeAvg11, UPrimeAvg22, UPrimeAvg33, UPrimeAvg12, UPrimeAvg13, UPrimeAvg23);

    // Set strain rate tensor using Uavg (or U) //
    volSymmTensorField Sij = symm(fvc::grad(U_));
    volSymmTensorField Sij_a = symm(fvc::grad(Uavg));

    // Get components of Sij //

    volScalarField S11 = Sij_a.component(tensor::XX);
    volScalarField S12 = Sij_a.component(tensor::XY);
    volScalarField S13 = Sij_a.component(tensor::XZ);
    volScalarField S22 = Sij_a.component(tensor::YY);
    volScalarField S23 = Sij_a.component(tensor::YZ);
    volScalarField S33 = Sij_a.component(tensor::ZZ);


    //Info << "\nDone composing strain rate tensor.\n" << endl;

    // Resolved turbulent production //

    RTP = UPrimeAvg11*S11 + UPrimeAvg12*S12 + UPrimeAvg13*S13 +
          UPrimeAvg22*S22 + UPrimeAvg23*S23 + UPrimeAvg33*S33;

    //Info << "\nDone calculating resolved turbulent production.\n" << endl;

    // Define Smagorinsky TKE //
    // The delta() method is user-specified
    ks = cI(D)*sqr(delta())*magSqr(D);

    //updateSubGridScaleFields(dev(Sij));

    // LES turbulent viscosity
    nuSgsL = cD(D)*sqr(delta())*sqrt(magSqr(D));


    // Define RANS viscosity (based on k-w SST) //
    //nuSgsR = a1_*k_/max(a1_*omega_,F2()*sqrt(2.0)*mag(symm(fvc::grad(Uavg))));
    Info<< "nuSgsR : min: " << min(nuSgsR) << " max: " << max(nuSgsR) << endl;

    volSymmTensorField tauSgs = nuSgsL*twoSymm(Sij);
    //Info << "\nDone composing model stresses.\n" << endl;

if (time > startAvg )
{
    if ( (time - dt) < startAvg )
    {
        tauSgsAvg = tauSgs;
    }
    else
    {
        tauSgsAvg = ( tauSgsAvg.oldTime()*(time-startAvg) + tauSgs*dt ) / (time-startAvg +dt);
         //Info << "\nDone averaging tauSgs.\n" << endl;
    }
}

    // tauRANS components //
    volScalarField tR11 = tauRANS.component(tensor::XX);
    volScalarField tR12 = tauRANS.component(tensor::XY);
    volScalarField tR13 = tauRANS.component(tensor::XZ);
    volScalarField tR22 = tauRANS.component(tensor::YY);
    volScalarField tR23 = tauRANS.component(tensor::YZ);
    volScalarField tR33 = tauRANS.component(tensor::ZZ);

    // tauSgs components //
    volScalarField tS11 = tauSgsAvg.component(tensor::XX);
    volScalarField tS12 = tauSgsAvg.component(tensor::XY);
    volScalarField tS13 = tauSgsAvg.component(tensor::XZ);
    volScalarField tS22 = tauSgsAvg.component(tensor::YY);
    volScalarField tS23 = tauSgsAvg.component(tensor::YZ);
    volScalarField tS33 = tauSgsAvg.component(tensor::ZZ);


    // Mean SGS Production //
    SP = tauSgsAvg && Sij_a;
    //Info << "\nDone computing Mean SGS production.\n" << endl;

    // Redefine flux to use Uavg instead of U
    surfaceScalarField phia_ = linearInterpolate(Uavg) & mesh_.Sf();

    const volScalarField S2(2*magSqr(symm(fvc::grad(Uavg))));
    volScalarField G("LESModel::G", nuSgsR*S2);

    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));


    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phia_, omega_)
      - fvm::Sp(fvc::div(phia_), omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)*S2
        //gamma(F1)*omega_/k_*nuSgsR*S2
        //min(gamma(F1)*S2, c1_*betaStar_*omega_*omega_)
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);
    Info<< "omega : min: " << min(omega_) << " max: " << max(omega_) << endl;

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phia_, k_)
      - fvm::Sp(fvc::div(phia_), k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G, c1_*betaStar_*k_*omega_)//Pk_tilde
      - fvm::Sp(betaStar_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);
    Info<< "k : min: " << min(k_) << " max: " << max(k_) << endl;

    nuSgsR = a1_*k_/max(a1_*omega_,F2()*sqrt(2.0)*mag(symm(fvc::grad(Uavg))));
    nuSgsR.correctBoundaryConditions();

    Info<< "nuSgsR : min: " << min(nuSgsR) << " max: " << max(nuSgsR) << endl;
     Info<< "nuSgsL : min: " << min(nuSgsL) << " max: " << max(nuSgsL) << endl;

    tauRANS = ((2.0/3.0)*I*k_) - nuSgsR*twoSymm(Sij_a); //already declared as a volSymmTensorField
    // RANS Production //
    RP = tauRANS && Sij_a;

    dimensionedScalar SMALLd = dimensionedScalar("1.0e-15", sqr(dimLength)/pow(dimTime,3), 1.0e-15);
    dimensionedScalar ONEd =  dimensionedScalar("1.0", sqr(dimLength)/pow(dimTime,3), 1.0);

    alpha = 1.0 - max(min((RTP/max(RP-SP, SMALLd)),1.0),1.0e-15);

    //Info << "\nDone computing alpha.\n" << endl;
    Info<< "alpha : min: " << gMin(alpha) << " max: " << gMax(alpha) << " average: " << gAverage(alpha) << endl;

    // Re-calculate viscosity
    // Using linear blending function
    // alpha = 1 - full RANS, 0 - full LES
    nuSgs_ = (1.0 - alpha)*nuSgsL + alpha*nuSgsR;
    nuSgs_.correctBoundaryConditions();
    //Info << "\nDone applying model blending.\n" << endl;
    ///////////////////////////////////

}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

