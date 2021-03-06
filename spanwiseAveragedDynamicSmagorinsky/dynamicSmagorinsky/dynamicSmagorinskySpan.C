/*---------------------------------------------------------------------------*\
dynamicSmagorinskySpan - Implementation of the spanwise averaged dynamic Smagorinsky
		     SGS model.
    
Copyright Information
    Local averaged dynamic smagorinsky model written by Dr. Alberto Passalacqua, Iowa State University 
    Modified by Zifei Yin, Iowa State University, for spanwise averaging of model constant, 2015
    zifeiyin@iastate.edu

License
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

#include "dynamicSmagorinskySpan.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicSmagorinskySpan, 0);
addToRunTimeSelectionTable(LESModel, dynamicSmagorinskySpan, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dynamicSmagorinskySpan::updateSubGridScaleFields
(
    const volSymmTensorField& D
)
{
    // The SGS viscosity is bounded so that nuEff cannot become negative.
    // Values are limited here, and not in nuEff, for consistency in stored
    // data and in submodels using nuSgs().
    // No warning message is printed when this limitation is applied.
    volScalarField averagedC_ = cD(D);
    averagedC_.internalField() =  meshIndexing1.collapse(cD(D));

    nuSgs_ = max(averagedC_*sqr(delta())*sqrt(magSqr(D)), -nu());
    nuSgs_.correctBoundaryConditions();
}

volScalarField dynamicSmagorinskySpan::cD
(
    const volSymmTensorField& D
) const
{
    tmp<volSymmTensorField> LL = 
	dev(filter_(sqr(U())) - (sqr(filter_(U()))));

    const volSymmTensorField MM
    (
        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D))
    );

    // Locally averaging MMMM on cell faces
    volScalarField MMMM = fvc::average(magSqr(MM));

    MMMM.max(VSMALL);

    // Performing local average on cell faces on return
    return 0.5*fvc::average(LL && MM)/MMMM;
}


volScalarField dynamicSmagorinskySpan::cI
(
    const volSymmTensorField& D
) const
{
    tmp<volScalarField> KK = 
	0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

    const volScalarField mm
    (
        sqr(delta())*(4*sqr(mag(filter_(D))) - filter_(sqr(mag(D))))
    );

    // Locally averaging mmmm on cell faces
    volScalarField mmmm = fvc::average(magSqr(mm));

    mmmm.max(VSMALL);

    // Performing local average on cell faces on return
    return fvc::average(KK*mm)/mmmm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicSmagorinskySpan::dynamicSmagorinskySpan
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(typeName, U, phi, transport),
    GenEddyVisc(U, phi, transport),

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
    spanwiseAverageDict
    (
        IOobject
        (
            "LESProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    meshIndexing1(mesh_, spanwiseAverageDict),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    updateSubGridScaleFields(dev(symm(fvc::grad(U))));
    

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void dynamicSmagorinskySpan::correct
(
    const tmp<volTensorField>& gradU
)
{
    LESModel::correct(gradU);

    const volSymmTensorField D(dev(symm(gradU)));

    k_ = cI(D)*sqr(delta())*magSqr(D);
    bound(k_,  kMin_);

    updateSubGridScaleFields(D);


    volScalarField kspan_ = k_;

    kspan_.internalField() =  meshIndexing1.collapse(k_);

}



bool dynamicSmagorinskySpan::read()
{
    if (GenEddyVisc::read())
    {
        filter_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
