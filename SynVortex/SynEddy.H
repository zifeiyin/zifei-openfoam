/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    SynEddyFvPatchVectorField

SourceFiles
    SynEddyFvPatchVectorField.C

Usage
    Synthetic eddy for LES inflow. 
    implemented according to Jarrin et al. International journal of heat
    and fluid flow, 2006

Author
    Zifei Yin, Iowa State University,  All rights reserved
    Email: zifeiyin@iastate.edu,  zifeiyin@gmail.com

sample of inflow boundary condition of U
    INFLOW
    {
        type                    SynEddy;
        freeStreamU             (1 0 0);
        refField                uniform (1 0 0);
        BoundIntensity          0.05;
        value                   uniform (1 0 0);
    }



\*---------------------------------------------------------------------------*/

#ifndef SynEddy_H
#define SynEddy_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class SynEddyFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class SynEddyFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Flow direction
        vector freeStreamU_;

        Field<vector>          refField_;

        Field<tensor>          Rstress_;


        scalar BoundIntensity_;






public:

    //- Runtime type information
    TypeName("SynEddy");


    // Constructors

        //- Construct from patch and internal field
        SynEddyFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SynEddyFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given SynEddyFvPatchVectorField
        //  onto a new patch
        SynEddyFvPatchVectorField
        (
            const SynEddyFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new SynEddyFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        SynEddyFvPatchVectorField
        (
            const SynEddyFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new SynEddyFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
