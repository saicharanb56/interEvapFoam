#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

immiscibleIncompressibleThreePhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
volScalarField& alpha3(mixture.alpha3());

const dimensionedScalar& rho1 = mixture.rho1(); //liquid
const dimensionedScalar& rho2 = mixture.rho2(); //vapour
const dimensionedScalar& rho3 = mixture.rho3(); //air

dimensionedScalar D23("D23", dimViscosity, mixture);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();

// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho1*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);


	// Create the diffusion coefficients for alpha2<->alpha3
        volScalarField Dc23(D23*max(alpha3, scalar(0))*pos0(alpha2));
        volScalarField Dc32(D23*max(alpha2, scalar(0))*pos0(alpha3));
        
        //const dimensionedScalar& recRho1 = 1/rho1;
        //const dimensionedScalar& recRho2 = 1/rho2;
        
        volScalarField limitedAlpha1(min(alpha1, scalar(1)));
        volScalarField limitedAlpha2(max(alpha2, scalar(0)));
        volScalarField limitedAlpha3(max(alpha3, scalar(0)));
        
        //const dimensionedScalar tempVar("1", dimensionSet(0, -1, 0, 0, 0), 1);
       
        const dimensionedScalar VSMALL("0.001", dimensionSet(1, -3, 0, 0, 0), 0.001);
              
volScalarField rhogp
(
    IOobject
    (
	"rhogp",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    (limitedAlpha2*rho2 + limitedAlpha3*rho3)/(limitedAlpha2 + limitedAlpha3 + scalar(0.001))
); 
        
volScalarField Y
(
    IOobject
    (
	"Y",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    (limitedAlpha2*rho2)/((rhogp + VSMALL)*(scalar(1.001) - limitedAlpha1))
);        

const dimensionedScalar tempValue("0.0001", dimless, 0.0001);

volVectorField mDot
(
    IOobject
    (
	"mDot",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
   (-tempValue*(Dc23 + Dc32)*rhogp*mag(fvc::grad(Y))*fvc::grad(limitedAlpha1)/(scalar(1) - Y)) 
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha2.name());

#include "createMRF.H"
