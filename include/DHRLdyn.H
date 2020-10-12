/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

Class
    Foam::incompressible::LESModels::DHRLdyn

Description
    The dynamic framework for hybrid RANS-LES modelling proposed by S. Bhushan
    and D. K. Walters (2011). 

    Herein, the RANS model k-w SST (Menter, 2003) is blended with the dynamic
    Smagorinsky LES model (Germano, 1991).

    The alpha coefficient serves as a blending function between RANS and LES
    stresses.

    alpha = 0 for full LES
    alpha = 1 for full RANS

SourceFiles
    DHRLdyn.C

\*---------------------------------------------------------------------------*/

#ifndef DHRLdyn_H
#define DHRLdyn_H

//#include "SSTEric/SSTEric.H"
//#include "Smagorinsky.H"
#include "LESModel.H"
#include "LESfilter.H"
//#include "LESModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class DHRLdyn Declaration
\*---------------------------------------------------------------------------*/

class DHRLdyn
:
    public LESModel
{
////
/*
           //- Update sub-grid scale fields
	   void updateSubGridScaleFields();
*/

protected:

    // Protected data

        // Model coefficients
            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar c1_;

            

	    //dimensionedScalar omegaMin_;//-- ERIC

	    dimensionedScalar ck_;
            dimensionedScalar ce_;
            
            
   
        //- Wall distance field
        //  Note: different to wall distance in parent LESModel
        wallDist y_;

        // Fields

            volScalarField k_;
            volScalarField ks;
            dimensionedScalar omegaMin_;
	    
            volScalarField omega_;
            volVectorField Uavg;
            volSymmTensorField tauSgsAvg;
            volSymmTensorField tauRANS;
            volScalarField RTP;
            volScalarField SP;
            volScalarField RP;
            volScalarField alpha;
            volScalarField nuSgs_;
            volScalarField nuSgsR;
            volScalarField nuSgsL;
            volScalarField UPrimeAvg11;
            volScalarField UPrimeAvg22;
            volScalarField UPrimeAvg33;
            volScalarField UPrimeAvg12;
            volScalarField UPrimeAvg13;
            volScalarField UPrimeAvg23;
            
            autoPtr<LESfilter> filterPtr_;
            LESfilter& filter_;
	            
	    	    


    // Protected Member Functions

        tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        tmp<volScalarField> F2() const;
        

        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaK(const volScalarField& F1) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        tmp<volScalarField> alphaOmega(const volScalarField& F1) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        tmp<volScalarField> beta(const volScalarField& F1) const
        {
            return blend(F1, beta1_, beta2_);
        }

        tmp<volScalarField> gamma(const volScalarField& F1) const
        {
            return blend(F1, gamma1_, gamma2_);
        }

        
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        //DHRLdyn(const DHRLdyn&);
        //DHRLdyn& operator=(const DHRLdyn&);

//- Update sub-grid scale fields
        void updateSubGridScaleFields(const volSymmTensorField& D);

        //- Calculate coefficients cD, cI from filtering velocity field
        dimensionedScalar cD(const volSymmTensorField& D) const;
        dimensionedScalar cI(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        DHRLdyn(const DHRLdyn&);
        DHRLdyn& operator=(const DHRLdyn&);


public:

    //- Runtime type information
    TypeName("DHRLdyn");


    // Constructors

        //- Construct from components
        DHRLdyn
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~DHRLdyn()
    {}


    // Member Functions

        //- Return the SGS turbulent kinetic energy.
        //virtual tmp<volScalarField> k() const;
          //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return omega
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }
/*
        virtual tmp<volScalarField> delta() const
        {
            return delta_;
        } */
          //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        /*
        tmp<volScalarField> ks(const tmp<volTensorField>& gradU) const
        {
            return (2.0*ck_/ce_)*sqr(delta())*magSqr(dev(symm(gradU)));
        } */
        //- Return the SGS turbulent disipation rate.
         virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

 //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                 new volScalarField("DkEff", alphaK(F1)*nuSgsR + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaOmega(F1)*nuSgsR + nu())
            );
        }
         

          //- Return the turb viscosity.
        virtual tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }


        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        virtual tmp<volSymmTensorField> devBeff() const;

        virtual tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
       

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

