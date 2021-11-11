{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    
    surfaceScalarField phir("phir", phic*mixture.nHatf());
    
    tmp<surfaceScalarField> talphaPhi1;
    tmp<surfaceScalarField> talphaPhi2;
    
    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
          - fvm::Sp(divU, alpha1) 
          ==
          -mag(mDot)/rho1
        );

        alpha1Eqn.solve();
        
        fvScalarMatrix alpha2Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha2)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha2)
          - fvm::Sp(divU, alpha2)
          - fvm::laplacian(Dc23 + Dc32, alpha2)
          ==
            mag(mDot)/rho2
        );

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
            
        Info<< "Phase-2 volume fraction = "
            << alpha2.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
            << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
            << endl;

        talphaPhi1 = alpha1Eqn.flux();
        talphaPhi2 = alpha2Eqn.flux();
    }
    
    volScalarField alpha10("alpha10", alpha1);
    volScalarField alpha20("alpha20", alpha2);
    
    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
       
        // Create the complete convection flux for alpha1
        tmp<surfaceScalarField> talphaPhiCorr1
        (                       
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha3, alpharScheme),
                alpha1,
                alpharScheme
            )
        );
        
         if (MULESCorr)
        {
            talphaPhiCorr1.ref() -= talphaPhi1();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talphaPhi1(),
                talphaPhiCorr1.ref(),
                (-mag(mDot))(),
               // (
               //     divU*(alpha10 - alpha100)
               // )(),
               zeroField(),
                1,
                0
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talphaPhi1.ref() += talphaPhiCorr1();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                talphaPhi1.ref() += 0.5*talphaPhiCorr1();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                talphaPhiCorr1.ref(),
                (-mag(mDot)/rho1)(),
                //(divU*alpha1)(),
                zeroField(),
                1,
                0
            );

            talphaPhi1 = talphaPhiCorr1;
        }

        // Create the complete flux for alpha2
        tmp<surfaceScalarField> talphaPhiCorr2
        (
            fvc::flux
            (
                phi,
                alpha2,
                alphaScheme
            )
            
          + fvc::flux
            (
                -fvc::flux(phir, alpha1, alpharScheme),
                alpha2,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            talphaPhiCorr2.ref() -= talphaPhi2();

            volScalarField alpha200("alpha200", alpha20);
            alpha20 = alpha2;

            MULES::correct
            (
                geometricOneField(),
                alpha2,
                talphaPhi2(),
                talphaPhiCorr2.ref(),
                (mag(mDot))(),
               // (
               //     divU*(alpha20 - alpha200)
               // )(),
               zeroField(),
                1,
                0
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talphaPhi2.ref() += talphaPhiCorr2();
            }
            else
            {
                alpha2 = 0.5*alpha2 + 0.5*alpha20;
                talphaPhi2.ref() += 0.5*talphaPhiCorr2();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha2,
                phi,
                talphaPhiCorr2.ref(),
                (mag(mDot)/rho2)(),
               // (divU*alpha2)(),
                zeroField(),
                1,
                0
            );

            talphaPhi2 = talphaPhiCorr2;
        }

        // Add the diffusive flux for alpha3->alpha2
        //alphaPhi2 -= fvc::interpolate(Dc32)*mesh.magSf()*fvc::snGrad(alpha1);
                 
        alpha3 = 1.0 - alpha1 - alpha2;
    }
    
    // Construct the complete mass flux
        rhoPhi =
              talphaPhi1()*(rho1 - rho3)
            + talphaPhi2()*(rho2 - rho3)
            + phi*rho3;

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
       // << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    Info<< "Vapour phase volume fraction = "
        << alpha2.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
       // << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
        << endl;
}
