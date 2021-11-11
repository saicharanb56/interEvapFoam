    MRF.correctBoundaryVelocity(U);
        
    volVectorField U_corrf //final U_corr defined at cell centre
    (
    	IOobject
    	(
	    "U_corrf",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
    	),
    	fvc::reconstruct(U_corr & mesh.Sf()) 
    );
    
    Info << "Min U_corrf = " << min(U_corrf) << endl;
    Info << "Max U_corrf = " << max(U_corrf) << endl;
    
    U -= U_corrf;
  
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::div(rhoPhi), U)
     // + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
     // + fvm::ddt(rho, U_corrf)
     // + fvm::div(rhoPhi, U_corrf)
     // - fvm::Sp(fvc::div(rhoPhi), U_corrf)  
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
