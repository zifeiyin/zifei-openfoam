/*---------------------------------------------------------------------------*\
dynamicSmagorinskySpan - Implementation of the dynamic Smagorinsky SGS model
		     as proposed by Lilly (1992) for OpenFOAM

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

Class
    Foam::incompressible::LESModels::dynamicSmagorinskySpan

Description
    The isochoric dynamic Smagorinsky model for incompressible flows.

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevails.
    Thus,
    @verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)

    where

        k = cI*delta^2*||D||^2
        nuSgs = cD*delta^2*||D||
        nuEff = nuSgs + nu

    In the dynamic version of the choric  Smagorinsky model
    the coefficients cI and cD are calculated during the simulation,

        cI=<K*m>_face/<m*m>_face

    and

        cD=<L.M>_face/<M.M>_face,

    where

        K = 0.5*(F(U.U) - F(U).F(U))
        m = delta^2*(4*||F(D)||^2 - F(||D||^2))
        L = dev(F(U*U) - F(U)*F(U))
        M = delta^2*(F(||D||*dev(D)) - 4*||F(D)||*F(dev(D)))
        <...>_face = face average
    @endverbatim

SourceFiles
    dynamicSmagorinskySpan.C
    
Authors
    Alberto Passalacqua <albertop@iastate.edu; albert.passalacqua@gmail.com>
    spanwise averaging is added by Zifei Yin  zifeiyin@gmail.com

References
    -	Lilly, D. K., A proposed modificaiton of the Germano subgrid-scale 
	closure method, Physics of Fluid A, 4 (3), 1992.

Notes
    Implementation of the dynamic Smagorinsky model with coefficients cD and
    cI computed as local average of their face values to avoid numerical 
    instabilities. 

    Negative values of the effective viscosity are removed by clipping it to
    zero (nuSgs is clipped to -nu)

    The code is known to work with OpenFOAM 2.3.0 and 2.3.x.

\*---------------------------------------------------------------------------*/

#ifndef dynamicSmagorinskySpan_H
#define dynamicSmagorinskySpan_H
#include "Average.H"
#include "Smagorinsky.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class dynamicSmagorinskySpan Declaration
\*---------------------------------------------------------------------------*/

class dynamicSmagorinskySpan
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;
        IOdictionary spanwiseAverageDict;
        meshIndex meshIndexing1;

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;
        
    //const dictionary& spanwiseAverageDict = U_.db().lookupObject<IOdictionary>
    //(
    // "LESProperties"
    //);
    //meshIndex meshIndexing1(mesh_, spanwiseAverageDict);

    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volSymmTensorField& D);
	
        //- Calculate coefficients cD, cI from filtering velocity field
        volScalarField cD(const volSymmTensorField& D) const;
        volScalarField cI(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        dynamicSmagorinskySpan(const dynamicSmagorinskySpan&);
        dynamicSmagorinskySpan& operator=(const dynamicSmagorinskySpan&);


public:

    //- Runtime type information
    TypeName("dynamicSmagorinskySpan");


    // Constructors

        //- Construct from components
        dynamicSmagorinskySpan
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );

    //- Destructor
    virtual ~dynamicSmagorinskySpan()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);


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
