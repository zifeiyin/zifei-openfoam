/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "dynamickOmegaDDES.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamickOmegaDDES, 0);
addToRunTimeSelectionTable(RASModel, dynamickOmegaDDES, dictionary);

volScalarField dynamickOmegaDDES::cD
(
    const volSymmTensorField& D,
    const volScalarField& w
) const
{
    
        
    tmp<volSymmTensorField> LL =
        dev(filter_(sqr(U())) - (sqr(filter_(U()))));

    const volSymmTensorField MM
    (
//        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D))
        sqr(delta())*(filter_(w*(D)) - 4*(filter_(w))*filter_(D))
    );

    // Locally averaging MMMM on cell faces
        volScalarField MMMM = fvc::average(magSqr(MM));
    
        MMMM.max(VSMALL);
    
        // Performing local average on cell faces on return
        return 0.5*fvc::average(LL && MM)/MMMM;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamickOmegaDDES::dynamickOmegaDDES
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
//            "betaStar",
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
//            "beta",
            "C_omega2",
            coeffDict_,
            0.072
        )
    ),
    alpha_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
//            "alpha",
            "C_omega1",
            coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            coeffDict_,
            0.5
        )
    ),

    Cu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cu",
            coeffDict_,
            0.12
         )
    ),
    Cd1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd1",
            coeffDict_,
            8.0
        )
    ),
    Cd2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd2",
            coeffDict_,
            3.0
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    betaV_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaV",
            coeffDict_,
            0.055
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),
    lDDES_
    (
        IOobject
        (
            "lDDES",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
//	dimensionedScalar("lDDES", dimLength*dimLength, 0.0)
	dimensionedScalar("lDDES", dimLength, 0.0)
    ),
    lRANS_
    (
        IOobject
        (
            "lRANS",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
//	dimensionedScalar("lRANS", dimLength*dimLength, 0.0)
	dimensionedScalar("lRANS", dimLength, 0.0)
    ),
    lLES_
    (
        IOobject
        (
            "lLES",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
//	dimensionedScalar("lLES", dimLength*dimLength, 0.0)
	dimensionedScalar("lLES", dimLength, 0.0)
    ),
    lzero_
    (
        IOobject
        (
            "lzero",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
//	dimensionedScalar("lzero", dimLength*dimLength, 0.0)
	dimensionedScalar("lzero", dimLength, 0.0)
    ),
    fd_
    (
        IOobject
        (
            "fd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
	dimensionedScalar("fd", dimVelocity/dimVelocity, 0.0)
    ),
    rd_
    (
        IOobject
        (
            "rd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
	dimensionedScalar("rd", dimVelocity/dimVelocity, 0.0)
    ),
    delta_
    (
        IOobject
        (
            "delta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
	dimensionedScalar("delta", dimLength, 1.0)
    ),
    cube_root_vol_
    (
        IOobject
        (
            "cube_root_vol",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
	dimensionedScalar("cube_root_vol", dimLength, 1.0)
    ),
    hmax_
    (
        IOobject
        (
            "hmax",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
	dimensionedScalar("hmax", dimLength, 1.0)
    ),

    CDES_
    (
        IOobject
        (
            "CDES",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
	dimensionedScalar("CDES", dimLength/dimLength, 0.12)
    ),
    CDES0_
    (
        IOobject
        (
            "CDES0",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
	dimensionedScalar("CDES0", dimLength/dimLength, 0.12)
    ),
        Lk_
    (
        IOobject
        (
            "Lk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
    mesh_,
    dimensionedScalar("Lk", dimLength, 1.0)
    ),
    Fr_
    (
        IOobject
        (
            "Fr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
    mesh_,
    dimensionedScalar("Fr", dimLength/dimLength, 1.0)
    ),
    Cbound_
    (
        IOobject
        (
            "Cbound",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_,
    dimensionedScalar("Cbound", dimLength/dimLength, 0.0)
    ),
    y_(mesh_),
    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    calcdelta();
    update_nut();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> dynamickOmegaDDES::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> dynamickOmegaDDES::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> dynamickOmegaDDES::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}

tmp<fvVectorMatrix> dynamickOmegaDDES::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

bool dynamickOmegaDDES::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        beta_.readIfPresent(coeffDict());
        alphaK_.readIfPresent(coeffDict());
        alphaOmega_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void dynamickOmegaDDES::update_nut()
{

    delta_ = fd_*cube_root_vol_ + (scalar(1) - fd_)*hmax_;

    volScalarField epsilon_tmp = sqr(CDES0_*hmax_)*omega_*scalar(2)*magSqr(symm(fvc::grad(U_)))  + Cmu_ * k_ * omega_;

    Lk_ = pow(nu(),scalar(3.0/4.0))/(pow(epsilon_tmp,scalar(1.0/4.0)));

    Fr_ = exp(-scalar(0.05)*hmax_/Lk_);

    Cbound_ = scalar(0.12) * (scalar(1.0)-tanh(Fr_*scalar(25.0)));

    lRANS_ = sqrt(k_)/omega_;

    CDES_ = sqrt(max(sqr(Cbound_ ), cD(dev(symm(fvc::grad(U_))), omega_ ) ) );

    lLES_ = CDES_ * delta_;
    lDDES_ = lRANS_ - fd_ * max(lzero_, lRANS_ - lLES_);
    nut_ = lDDES_ * lDDES_  * omega_;

    nut_.correctBoundaryConditions();
}

void dynamickOmegaDDES::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    const volScalarField O2(2*magSqr(skew(fvc::grad(U_))));
//    rd_ = (nut_ + nu())/(sqr(kappa_) * sqr(y_) * sqrt(scalar(0.5)*(S2 + O2)));
//--------------k/omega instead of nut for rd---------------------------
    rd_ = ((k_/omega_) + nu())/(sqr(kappa_) * sqr(y_) * sqrt(scalar(0.5)*(S2 + O2)));
    fd_ = scalar(1) - tanh(pow(Cd1_*rd_, Cd2_));
    // Re-calculate viscosity
    update_nut();

    volScalarField G("dynamickOmegaDDES:G", nut_*2*magSqr(symm(fvc::grad(U_))));

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(), omega_)
     ==
//        alpha_*G*omega_/k_
        alpha_*2*magSqr(symm(fvc::grad(U_)))
      - fvm::Sp(beta_*omega_, omega_)
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);
    //need the next line uncomment to prevent negative omega if bad initial condition is used
	//omega_ = max(omega_, scalar(0.02)*sqrt(2*magSqr(symm(fvc::grad(U_))))); // for preventing negative omega to ensure stability. Doesn't affect simulation result.


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(Cmu_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);

    Info << "min CDES = " << min(CDES_)<< "  max CDES = " << max(CDES_) << endl;

}

void dynamickOmegaDDES::calcdelta()
{
	label nD = mesh_.nGeometricD();
    if(nD==2)
    {
        WarningIn("dynamickOmegaDDES::calcDelta()")
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;
    }

    const cellList& cells = mesh_.cells();
    forAll(cells,cellI)
    {
        scalar deltaMaxTmp = 0.0;
  //      scalar cellLengthTmp = 0.0;
        //scalar faceNumberCount = 0.0;
        const labelList& cFaces = mesh_.cells()[cellI];
        const point& centrevector = mesh_.cellCentres()[cellI];

        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& facevector = mesh_.faceCentres()[faceI];
            scalar tmp = 2*mag(facevector - centrevector);
            deltaMaxTmp = max(deltaMaxTmp, tmp);

    //        cellLengthTmp = cellLengthTmp + magSqr(tmp);
            //faceNumberCount = faceNumberCount + scalar(1.0);

        }

        hmax_[cellI] = deltaMaxTmp;
//        Lc_[cellI] = sqrt(cellLengthTmp/scalar(2.0));
    }

    cube_root_vol_.internalField() = pow(mesh_.V(), 1.0/3.0);


}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
