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

\*---------------------------------------------------------------------------*/
#include "SynEddy.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "Random.H"
//#include "fvMesh.H"
//#include "fvCFD.H"
//#include "surfaceFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Foam::Random ranGen(Foam::label(time(0)));

Foam::scalar xmax_ = -1e6;
Foam::scalar xmin_ = 1e6;
Foam::scalar ymax_ = -1e6;
Foam::scalar ymin_ = 1e6;
Foam::scalar zmax_ = -1e6;
Foam::scalar zmin_ = 1e6;
Foam::scalar margin_ = 0.005;

//Foam::vector freeStreamU_ = Foam::vector(1.0,0.0,0.0);

const int eddyNum_ = 100;
int kcount_ = 0;
Foam::vector eddyCoord_[eddyNum_];
//Foam::scalar eddyIntensity_[eddyNum_];
Foam::vector randomSign_[eddyNum_];

bool setBox = false;
bool eddyGenerated = false;

Foam::scalar eddyConvectionTime_ = 0.0;

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SynEddyFvPatchVectorField::SynEddyFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    freeStreamU_(1, 0, 0),
    refField_(p.size()),
    Rstress_(),
    BoundIntensity_(5e-4)
{}


SynEddyFvPatchVectorField::SynEddyFvPatchVectorField
(
    const SynEddyFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    freeStreamU_(ptf.freeStreamU_),
    refField_(ptf.refField_, mapper),
    Rstress_(ptf.Rstress_, mapper),
    BoundIntensity_(ptf.BoundIntensity_)
{}


SynEddyFvPatchVectorField::SynEddyFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    freeStreamU_(dict.lookup("freeStreamU")),
    refField_("refField", dict, p.size()),
    Rstress_(p.size(), pTraits<tensor>::zero),
    BoundIntensity_(readScalar(dict.lookup("BoundIntensity")))
{
    evaluate();
}


SynEddyFvPatchVectorField::SynEddyFvPatchVectorField
(
    const SynEddyFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    freeStreamU_(fcvpvf.freeStreamU_),
    refField_(fcvpvf.refField_),
    Rstress_(fcvpvf.Rstress_),
    BoundIntensity_(fcvpvf.BoundIntensity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SynEddyFvPatchVectorField::updateCoeffs()
{

    if (updated())
    {
        return;
    }

    //int oldTag = UPstream::msgType();
    //UPstream::msgType() = oldTag+1;

    //Field<vector>& patchField = *this;

    scalar t_ = this->db().time().timeOutputValue();
    scalar startT_ = this->db().time().startTime().value();
    scalar deltaT_ = db().time().deltaTValue();

if (Pstream::myProcNo()==0)
{
    if (t_ > startT_ + 1.5 * deltaT_)
    {
        const vectorField& c = patch().Cf();
    
        // setting synthetic eddy box
        // checking max and minimum x, y, and z components 
        // in patch and then apply margin
        if (setBox==false)
        {
            forAll(c,cellI)
            {
                xmin_ = min(xmin_,c[cellI].component(vector::X));
                xmax_ = max(xmax_,c[cellI].component(vector::X));
                ymin_ = min(ymin_,c[cellI].component(vector::Y));
                ymax_ = max(ymax_,c[cellI].component(vector::Y));
                zmin_ = min(ymin_,c[cellI].component(vector::Z));
                zmax_ = max(ymax_,c[cellI].component(vector::Z));
            }
            xmax_ += margin_;
            xmin_ -= margin_;
            ymax_ += margin_;
            ymin_ -= margin_;
            zmax_ += margin_;
            zmin_ -= margin_;
            setBox = true;
            Info << "setting bounding box for synthetic eddies" << endl;
            Info << "x = " << xmin_ << " to " << xmax_ << endl;
            Info << "y = " << ymin_ << " to " << ymax_ << endl;
            Info << "z = " << zmin_ << " to " << zmax_ << endl;
        }
    
        //random number generator to generate eddy coordinates and intensities
        if (eddyGenerated==false)
        {
            Info << "Generating random eddy coordinates" << endl;
            for (int i=0; i<eddyNum_; i++)
            {
                eddyCoord_[i].component(vector::X) = xmin_ + (xmax_ - xmin_) * ranGen.scalar01();
                eddyCoord_[i].component(vector::Y) = ymin_ + (ymax_ - ymin_) * ranGen.scalar01();
                eddyCoord_[i].component(vector::Z) = zmin_ + (zmax_ - zmin_) * ranGen.scalar01();
                randomSign_[i][0] = sign(ranGen.scalar01()-scalar(0.5));
                randomSign_[i][1] = sign(ranGen.scalar01()-scalar(0.5));
                randomSign_[i][2] = sign(ranGen.scalar01()-scalar(0.5));
            }
            eddyGenerated = true;
        }
    
        //coordinates of cell center
        //initialize velocity field
        const Field<vector> cellloc = patch().Cf();
        const fvMesh& mesh = patch().boundaryMesh().mesh();
    
        volVectorField Ufield = mesh.lookupObject<volVectorField>("U");
        //volVectorField Ufield = mesh.lookupObject<volVectorField>("U");
    
        Ufield.internalField() = refField_;//freeStreamU_;
        volVectorField Uturbulence = Ufield - Ufield; // make initial fluctuation to be 0

        // convect eddies in the domain
        if (eddyConvectionTime_ < t_)
        {
           //for (int i=0; i<eddyNum_; i++)
           //{
           //    Info << eddyCoord_[i] << endl;
           //}

            for (int i=0; i<eddyNum_; i++)
            {
                //eddyCoord_[i] += freeStreamU_ * deltaT_;
                eddyCoord_[i] += freeStreamU_ * deltaT_;
                kcount_ = 0; // counting updated eddy mumber
                //convect back to inflow plane when it is out of margin
                if (eddyCoord_[i].component(vector::X) > margin_)
                {
                    eddyCoord_[i].component(vector::X) -= scalar(2)*margin_;
                    eddyCoord_[i].component(vector::Y) = ymin_ + (ymax_ - ymin_) * ranGen.scalar01();
                    eddyCoord_[i].component(vector::Z) = zmin_ + (zmax_ - zmin_) * ranGen.scalar01();
                    randomSign_[i][0] = sign(ranGen.scalar01()-scalar(0.5));
                    randomSign_[i][1] = sign(ranGen.scalar01()-scalar(0.5));
                    randomSign_[i][2] = sign(ranGen.scalar01()-scalar(0.5));
                    kcount_ += 1;
                }
            }
            eddyConvectionTime_ += deltaT_;
            if (kcount_ > 0)
            {
                Info << kcount_ << "  eddies regenerated" << endl;
            }
            else
            {
                Info << "  eddies convected" << endl;
            }
        }
    
        // compute uprime vprime wprime
        //standard normal distribution
        //sigma = 1
        forAll(refField_,cellI)
        {
            for (int i = 0; i<eddyNum_; i++)
            {
                for (int j = 0; j < 3 ; j++)
                {
                    Uturbulence[cellI].component(j) += BoundIntensity_*mag( freeStreamU_ )
                        * scalar(1.0)/Foam::sqrt(scalar(eddyNum_))
                        *randomSign_[i][j]/(Foam::sqrt(scalar(6.2832))) 
                        *exp(-scalar(0.5)
                                *(
                                    magSqr(eddyCoord_[i] - cellloc[cellI])
                                    /(margin_)//normalization of x axis
                                )
                            );
                }
            }
        }
       
            //reconstruct Uturbulence using R11 R22 R33
        //for isotropic turbulence, there is no need.
        //Info << "perturbation magnitude =" << Uturbulence[0] << endl;
        vectorField::operator=(refField_+Uturbulence);
    } //t > 1.5deltaT + startT
    else
    {
        eddyConvectionTime_ = startT_+ deltaT_;
        vectorField::operator==(refField_);
        Info << "Using reference field" << endl;
    }

}//if (rank==0)


//UPstream::msgType() = oldTag;
//fixedValueFvPatchVectorField::updateCoeffs();
}//updateCoeffs()




// Write
void SynEddyFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("freeStreamU")
        << freeStreamU_ << token::END_STATEMENT << nl;
    refField_.writeEntry("refField", os);
    os.writeKeyword("BoundIntensity")
        << BoundIntensity_ << token::END_STATEMENT << nl;
 this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, SynEddyFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
