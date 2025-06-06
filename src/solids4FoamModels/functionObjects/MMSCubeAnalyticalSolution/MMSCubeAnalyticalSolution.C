/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*----------------------------------------------------------------------------*/

#include "MMSCubeAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MMSCubeAnalyticalSolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        MMSCubeAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::symmTensor Foam::MMSCubeAnalyticalSolution::MMSCubeStress
(
    const vector& point
)
{
    symmTensor sigma = symmTensor::zero;

    // Shear modulus
    const scalar mu = E_/(2.0*(1.0 + nu_));

    // Lambda parameter
    const scalar lambda = (E_*nu_)/((1.0+nu_)*(1.0-2.0*nu_));

    //pi
    const scalar pi = Foam::constant::mathematical::pi;

    sigma.xx() =
    lambda*(4*ax_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z())
    + 2*ay_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z())
    + az_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))
    + 8*ax_*mu*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    sigma.yy() =
	lambda*(4*ax_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z())
	+ 2*ay_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z())
	+ az_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))
	+ 4*ay_*mu*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z());

    sigma.zz() =
	lambda*(4*ax_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z())
	+ 2*ay_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z())
	+ az_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))
	+ 2*az_*mu*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y());

    sigma.xy() =
	2*ax_*mu*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z())
	+ 4*ay_*mu*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    sigma.yx() = sigma.xy();

    sigma.yz() =
    ay_*mu*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())
    + 2*az_*mu*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z());

    sigma.zy() = sigma.yz();

    sigma.xz() =
    ax_*mu*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())
    + 4*az_*mu*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    sigma.zx() = sigma.xz();

    return sigma;
}

Foam::symmTensor Foam::MMSCubeAnalyticalSolution::MMSCubeEpsilon
(
    const vector& point
)
{
    symmTensor epsilon = symmTensor::zero;

    //pi
    const scalar pi = Foam::constant::mathematical::pi;

    epsilon.xx() =
    4*ax_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    epsilon.yy() =
    2*ay_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z());

    epsilon.zz() =
    az_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y());

    epsilon.xy() =
    ax_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z()) + 2*ay_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    epsilon.yx() = epsilon.xy();

    epsilon.yz() =
    (ay_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))/2 + az_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z());

    epsilon.zy() = epsilon.yz();

    epsilon.xz() =
    (ax_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))/2 + 2*az_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    epsilon.zx() = epsilon.xz();

    return epsilon;
}

Foam::vector Foam::MMSCubeAnalyticalSolution::MMSCubeDisplacement
(
    const vector& point
)
{
    //pi
    const scalar pi = Foam::constant::mathematical::pi;

    return vector
    (
        ax_*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z()),
        ay_*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z()),
        az_*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z())
    );
}

bool Foam::MMSCubeAnalyticalSolution::writeData()
{
    // Lookup the solid mesh
    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    // Lookup the point mesh
    const pointMesh& pMesh = mesh.lookupObject<pointMesh>("pointMesh");

    // Point coordinates
    const pointField& points = mesh.points();

    //Cell-centre coordinates
    const volVectorField& C = mesh.C();
    const vectorField& CI = C;

    // Point analytical fields
    if (pointDisplacement_ || pointStress_ || pointEpsilon_)
    {
        volSymmTensorField analyticalStress
        (
            IOobject
            (
                "analyticalCellStress",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
            "calculated"
        );

        pointSymmTensorField analyticalPointEpsilon
        (
            IOobject
            (
                "analyticalPointEpsilon",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
            "calculated"
        );

		pointScalarField analyticalPointEpsilonEq
		(
		    IOobject
		    (
		        "analyticalPointEpsilonEq",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
		    ),
		    pMesh,
		    dimensionedScalar("0", dimless, 0.0),
		    "calculated"
		);

        volScalarField analyticalStressEq
        (
            IOobject
            (
                "analyticalCellStressEq",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPressure, 0),
            "calculated"
        );

        pointVectorField analyticalPointD
        (
            IOobject
            (
                "analyticalPointD",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("zero", dimLength, vector::zero)
        );

        volVectorField analyticalD
        (
            IOobject
            (
                "analyticalD",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimLength, vector::zero)
        );

        symmTensorField& sI = analyticalStress;
        symmTensorField& pEI = analyticalPointEpsilon;
        scalarField& pEEqI = analyticalPointEpsilonEq;
        scalarField& sEqI = analyticalStressEq;
        vectorField& aPDI = analyticalPointD;
        vectorField& aDI = analyticalD;


        forAll(sI, cellI)
        {
            if (pointStress_)
            {
                sI[cellI] = MMSCubeStress(CI[cellI]);
                sEqI[cellI] = sqrt((sqr(sI[cellI].xx() - sI[cellI].yy()) + sqr(sI[cellI].yy() - sI[cellI].zz()) + sqr(sI[cellI].zz() - sI[cellI].xx()) +
                6*(sqr(sI[cellI].xy()) + sqr(sI[cellI].yz()) + sqr(sI[cellI].zx())))/2);
            }

        }

        forAll(pEI, pointI)
        {
            if (pointEpsilon_)
            {
                pEI[pointI] = MMSCubeEpsilon(points[pointI]);
                pEEqI[pointI] = sqrt((2.0/3.0)*magSqr(dev(pEI[pointI])));

            }

        }

        forAll(aPDI, pointI)
        {

            if (pointDisplacement_)
            {
                aPDI[pointI] = MMSCubeDisplacement(points[pointI]);
            }
        }

        forAll(aDI, cellI)
        {

            if (pointDisplacement_)
            {
                aDI[cellI] = MMSCubeDisplacement(CI[cellI]);
            }
        }

        forAll(analyticalStress.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
#ifdef OPENFOAM_NOT_EXTEND
                symmTensorField& sP = analyticalStress.boundaryFieldRef()[patchI];
                scalarField& sEqP = analyticalStressEq.boundaryFieldRef()[patchI];
#else
                symmTensorField& sP = analyticalStress.boundaryField()[patchI];
                scalarField& sEqP = analyticalStressEq.boundaryField()[patchI];
#endif
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(sP, faceI)
                {
                	sP[faceI] = MMSCubeStress(CP[faceI]);
                	sEqP[faceI] = sqrt((sqr(sP[faceI].xx() - sP[faceI].yy()) + sqr(sP[faceI].yy() - sP[faceI].zz()) + sqr(sP[faceI].zz() - sP[faceI].xx()) +
                6*(sqr(sP[faceI].xy()) + sqr(sP[faceI].yz()) + sqr(sP[faceI].zx())))/2);
                }
            }
        }

        // Write point analytical fields
        if (pointStress_)
        {
            Info<< "Writing analyticalPointStress and analyticalPointStressEq"
                << nl << endl;
            analyticalStress.write();
            analyticalStressEq.write();
        }

        if (pointEpsilon_)
        {
            Info<< "Writing analyticalPointEpsilon and analyticalPointEpsilonEq"
                << nl << endl;
            analyticalPointEpsilon.write();
            analyticalPointEpsilonEq.write();
        }

        if (pointDisplacement_)
        {
            Info<< "Writing analyticalPointDisplacement and analyticalDisplacement"
                << nl << endl;
            analyticalPointD.write();
            analyticalD.write();
        }

        if
        (
            pointDisplacement_
         && mesh.foundObject<pointVectorField>("pointD")
        )
        {
            const pointVectorField& pointD =
                mesh.lookupObject<pointVectorField>("pointD");

            const pointVectorField diff
            (
                "pointDDifference", analyticalPointD - pointD
            );
            Info<< "Writing pointDDifference field" 
		<< nl << endl;
            diff.write();

            // const vectorField& diffI = diff;
            // Info<< "    Displacement norms: mean L1, mean L2, LInf: " << nl
            //     << "	Magnitude: " << gAverage(mag(diffI))
            //     << " " << Foam::sqrt(gAverage(magSqr(diffI)))
            //     << " " << gMax(mag(diffI))
            //     << nl
            //     << "	u: " << gAverage(mag(diffI.component(0)))
            //     << " " << Foam::sqrt(gAverage(magSqr(diffI.component(0))))
            //     << " " << gMax(mag(diffI.component(0)))
            //     << nl
            //     << "	v: " << gAverage(mag(diffI.component(1)))
            //     << " " << Foam::sqrt(gAverage(magSqr(diffI.component(1))))
            //     << " " << gMax(mag(diffI.component(1)))
            //     << nl
            //     << "	w: " << gAverage(mag(diffI.component(2)))
            //     << " " << Foam::sqrt(gAverage(magSqr(diffI.component(2))))
            //     << " " << gMax(mag(diffI.component(2)))
            //     << nl << endl;
        }

        if
        (
            pointDisplacement_
         && mesh.foundObject<volVectorField>("D")
        )
        {
            const volVectorField& D =
                mesh.lookupObject<volVectorField>("D");

            const volVectorField diff
            (
                "DDifference", analyticalD - D
            );
            Info<< "Writing DDifference field" 
		<< nl << endl;
            diff.write();

            // const vectorField& diffI = diff;
            //  Info<< "    Displacement norms: mean L1, mean L2, LInf: " << nl
            //     << "	Magnitude: " << gAverage(mag(diffI))
            //     << " " << Foam::sqrt(gAverage(magSqr(diffI)))
            //     << " " << gMax(mag(diffI))
            //     << nl
            //     << "	u: " << gAverage(mag(diffI.component(0)))
            //     << " " << Foam::sqrt(gAverage(magSqr(diffI.component(0))))
            //     << " " << gMax(mag(diffI.component(0)))
            //     << nl
            //     << "	v: " << gAverage(mag(diffI.component(1)))
            //     << " " << Foam::sqrt(gAverage(magSqr(diffI.component(1))))
            //     << " " << gMax(mag(diffI.component(1)))
            //     << nl
            //     << "	w: " << gAverage(mag(diffI.component(2)))
            //     << " " << Foam::sqrt(gAverage(magSqr(diffI.component(2))))
            //     << " " << gMax(mag(diffI.component(2)))
            //     << nl << endl;
        }

        if
        (
            pointEpsilon_
         && mesh.foundObject<pointSymmTensorField>("pEpsilon")
        )
        {
            const pointSymmTensorField& pointEpsilon =
                mesh.lookupObject<pointSymmTensorField>("pEpsilon");

            const pointSymmTensorField diffEpsilon
            (
                "pointEpsilonDifference", analyticalPointEpsilon - pointEpsilon
            );
            Info<< "Writing pointPointEpsilonDifference field" 
		<< nl << endl;
            diffEpsilon.write();

            // const symmTensorField& diffEpsilonI = diffEpsilon;
            // Info<< "    pEpsilon norms: mean L1, mean L2, LInf: " << nl
            //    << "    Component XX: " << gAverage(mag(diffEpsilonI.component(0)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(0))))
            //    << " " << gMax(mag(diffEpsilonI.component(0)))
            //    << nl
            //    << "    Component XY: " << gAverage(mag(diffEpsilonI.component(1)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(1))))
            //    << " " << gMax(mag(diffEpsilonI.component(1)))
            //    << nl
            //    << "    Component XZ: " << gAverage(mag(diffEpsilonI.component(2)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(2))))
            //    << " " << gMax(mag(diffEpsilonI.component(2)))
            //    << nl
            //    << "    Component YY: " << gAverage(mag(diffEpsilonI.component(3)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(3))))
            //    << " " << gMax(mag(diffEpsilonI.component(3)))
            //    << nl
            //    << "    Component YZ: " << gAverage(mag(diffEpsilonI.component(4)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(4))))
            //    << " " << gMax(mag(diffEpsilonI.component(4)))
            //    << nl
            //    << "    Component ZZ: " << gAverage(mag(diffEpsilonI.component(5)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(5))))
            //    << " " << gMax(mag(diffEpsilonI.component(5)))
            //    << nl << endl;
        }


        if
        (
            pointStress_
         && mesh.foundObject<volSymmTensorField>("sigma")
        )
        {

            const volSymmTensorField& pointSigma =
                mesh.lookupObject<volSymmTensorField>("sigma");

//            const scalarField& pointSigmaEq(mesh.nPoints(),0);

//            forAll (pointSigma, pointI)
//            {
//            pointSigmaEq[pointI] =
//            sqrt((sqr(pointSigma[pointI].xx() - pointSigma[pointI].yy()) + sqr(pointSigma[pointI].yy() - pointSigma[pointI].zz()) +
//            sqr(pointSigma[pointI].zz() - pointSigma[pointI].xx()) + 6*(sqr(pointSigma[pointI].xy()) + sqr(pointSigma[pointI].yz())
//            + sqr(pointSigma[pointI].zx())))/2);
//            }


            const volSymmTensorField diffSigma
            (
                "pointSigmaDifference", analyticalStress - pointSigma
            );

//            const volScalarField diffSigmaEq
//            (
//                "pointSigmaEqDifference", analyticalStressEq - pointSigmaEq
//            );

            Info<< "Writing pointSigmaDifference and pointSigmaEqDifference fields" 
		<< nl << endl;
            diffSigma.write();
            //diffSigmaEq.write();

            // const symmTensorField& diffSigmaI = diffSigma;
            // const scalarField& diffSigmaEqI = diffSigmaEq;
            // Info<< "    Sigma norms: mean L1, mean L2, LInf: " << nl
            //    << "    Component XX: " << gAverage(mag(diffSigmaI.component(0)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffSigmaI.component(0))))
            //    << " " << gMax(mag(diffSigmaI.component(0)))
            //    << nl
            //    << "    Component XY: " << gAverage(mag(diffSigmaI.component(1)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffSigmaI.component(1))))
            //    << " " << gMax(mag(diffSigmaI.component(1)))
            //    << nl
            //    << "    Component XZ: " << gAverage(mag(diffSigmaI.component(2)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffSigmaI.component(2))))
            //    << " " << gMax(mag(diffSigmaI.component(2)))
            //    << nl
            //    << "    Component YY: " << gAverage(mag(diffSigmaI.component(3)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffSigmaI.component(3))))
            //    << " " << gMax(mag(diffSigmaI.component(3)))
            //    << nl
            //    << "    Component YZ: " << gAverage(mag(diffSigmaI.component(4)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffSigmaI.component(4))))
            //    << " " << gMax(mag(diffSigmaI.component(4)))
            //    << nl
            //    << "    Component ZZ: " << gAverage(mag(diffSigmaI.component(5)))
            //    << " " << Foam::sqrt(gAverage(magSqr(diffSigmaI.component(5))))
            //    << " " << gMax(mag(diffSigmaI.component(5)))
            //    << nl << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MMSCubeAnalyticalSolution::MMSCubeAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu"))),
    pointDisplacement_
    (
        dict.lookupOrDefault<Switch>("pointDisplacement", true)
    ),
    pointStress_
    (
        dict.lookupOrDefault<Switch>("pointStress", true)
    ),
    pointEpsilon_
    (
        dict.lookupOrDefault<Switch>("pointEpsilon", true)
    ),
    ax_(readScalar(dict.lookup("ax"))),
    ay_(readScalar(dict.lookup("ay"))),
    az_(readScalar(dict.lookup("az")))
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (E_ < SMALL || nu_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "E and nu should be positive!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::MMSCubeAnalyticalSolution::start()
{
    return true;
}


#if FOAMEXTEND
    bool Foam::MMSCubeAnalyticalSolution::execute(const bool forceWrite)
#else
    bool Foam::MMSCubeAnalyticalSolution::execute()
#endif
{
    return writeData();
}


bool Foam::MMSCubeAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::MMSCubeAnalyticalSolution::write()
{
    return true;
}
#endif

// ************************************************************************* //
