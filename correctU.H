//surfaceVectorField U_rho(fvc::interpolate(U));

surfaceScalarField fA(fvc::interpolate(alpha1)); //portion of a face that is wetted

volScalarField temp //temporary field     
    (             
    IOobject 
        (                     
        "temp",
        mesh,
        IOobject::NO_READ,                     
        IOobject::AUTO_WRITE
    	), 
    	mesh,
    	dimensionedScalar("one",dimless,1.0) 
    );
    
const surfaceScalarField recRhof(fvc::interpolate(temp*((1/rho2)-(1/rho1))));
surfaceVectorField mDotf(fvc::interpolate(mDot)); //face value of mDot vector

label celli;
const cellList& Cells = mesh.cells(); //list of all cell IDs
const faceList& faces = mesh.faces(); //list of all face IDs

for(celli = 0; celli < mesh.nCells(); celli++)
{
    if(alpha1[celli] > 0 && alpha1[celli] < 1)
    {	
    	scalar pValue = p.ref()[celli];
    	
    	Y[celli] = temp[celli];
    	Y[celli] = Y[celli]*1000000*18;
    	Y[celli] = Y[celli]/(100000*18 + (pValue - 1000000)*29);
    }
}

surfaceScalarField rhof(fvc::interpolate(rho));
surfaceScalarField rho1f(fvc::interpolate(temp*rho1));
surfaceScalarField rhogpf(fvc::interpolate(rhogp));

surfaceVectorField Uri(fvc::interpolate(U));
surfaceVectorField Ugpi(fvc::interpolate(U));

surfaceScalarField eL //edge length field defined at a cell face   
    (             
    IOobject 
    	(                     
        "eL",
        mesh,
        IOobject::NO_READ,                     
        IOobject::AUTO_WRITE
        ), 
    	mesh,
    	dimensionedScalar("one",dimless, 0.001) 
    );

forAll(eL, facei)
{
	//dimensionedScalar edgeLength("edgeLength", dimLength, mesh.V()[mesh.faceOwner()[facei]]/mesh.magSf()[facei]);
	
	scalar edgeLength = mesh.V()[mesh.faceOwner()[facei]]/mesh.magSf()[facei];
	
	eL[facei] *= edgeLength;
	
	Uri[facei] -= mDotf[facei]*recRhof[facei]*eL[facei]*(rhogpf[facei]/rhof[facei]);
	
	Ugpi[facei] += mDotf[facei]*recRhof[facei]*eL[facei]*(rho1f[facei]/rhof[facei])*(fA[facei]/(1.001 - fA[facei]));
}

//STEP 2

surfaceVectorField U_vol = Uri*fA + Ugpi*(1-fA);

//surfaceVectorField U_vol = U_rho 
//			 + (mDotf*recRhof*eL*(fvc::interpolate(temp*rho1)/fvc::interpolate(rho)) 
//			 - mDotf*recRhof*eL*(fvc::interpolate(rhogp)/fvc::interpolate(temp*rho)))*fA; //volume weighted velocity


const dimensionedScalar VSMALL2("0.001", dimensionSet(1, -3, -1, 0, 0), 0.001);

volVectorField Err = mDot*((1/rho2)-(1/rho1)) - fvc::div(U_vol & mesh.Sf())*mDot/(mag(mDot) + VSMALL2); //Calculation of Error
	
//STEP 3

labelList exclude(mesh.nCells()); //list of cell IDs to exclude as those cells would have already been assigned credit scores

forAll(exclude, celli)  //initialising the exclude list
{
    exclude[celli] = -1;
}

label i = 0; 
scalar maxAlpha = 0; //max value in a particular cell
    
while(i < mesh.nCells())
{   
    label maxCelli = -1;
    
    for (celli = 0; celli < mesh.nCells(); celli++)
    {	
    	if(alpha1[celli]>0 && alpha1[celli]<1)
    	{
    	    	bool calc; //variable to ensure that maxCelli of last iter is excluded in current iter
    		label j = 0;
    	
    		while(j<mesh.nCells()) //loop to ensure that the maxCelli of the last iter is excluded in the current iter
    		{
    			if(celli == exclude[j])
    			{
    				calc = false;
    				break;
       			}
    		
    			else
    			{
    				calc = true;
    			}
    		
    			j++;
    		}
    	    				
    	  	if(calc)
    	  	{
    	  		if(alpha1[celli] > maxAlpha)
    			{
    				maxAlpha = alpha1[celli];
    				maxCelli = celli;
    			}
    		}
    	}
    	    
    }
               
    exclude[i] = maxCelli; //also the list that contains the cellIDs in descending order of their alpha1 i.e. in the order they need to be solved
    i++;
}

for(i = 0; i < mesh.nCells(); i++)
{   
    label facei = 0;
    
    if(exclude[i] == -1)
    {
    	break;
    }
    
    while(facei < mesh.nFaces())
    {
    	if(alpha1[mesh.faceOwner()[facei]] > 0 && alpha1[mesh.faceOwner()[facei]] < 1 )
    	{
    	  	if(mesh.faceOwner()[facei] == exclude[i])
    		{
    			celli = exclude[i];
    			const cell& Cell = Cells[celli];
    			
    			forAll(U_vol, facei)
    			{
    				label oppFacei = Cell.opposingFaceLabel(facei, faces);
    				
    				if(facei != -1 && oppFacei != -1)
    				{
    				U_vol[facei] += Err[celli]*eL[facei]*(1/rhof[facei])/((1/rhof[facei]) + (1/rhof[oppFacei]));
    										
    				U_vol[oppFacei] += Err[celli]*eL[oppFacei]*(1/rhof[oppFacei])/((1/rhof[facei]) + (1/rhof[oppFacei]));
    				}
    		   	}
    		}
    	}
    	
    	facei++;
    }
}

surfaceVectorField Ur(U_vol);
surfaceVectorField Ugp(U_vol);

forAll(eL, facei)
{
	scalar edgeLength = mesh.V()[mesh.faceOwner()[facei]]/mesh.magSf()[facei];
	
	eL[facei] *= edgeLength;
	
	Ur[facei] += mDotf[facei]*recRhof[facei]*eL[facei]*(fA[facei]/(1.001 - fA[facei])); //interface velocity defined at a cell face

	Ugp[facei] -= mDotf[facei]*recRhof[facei]*eL[facei]; //gas phase velocity defined at a cell face
}

//STEP 4
surfaceVectorField U_corr = (Ur*fA*rho1 + Ugp*(1-fA)*fvc::interpolate(rhogp))/fvc::interpolate(rho); //corrected mass weighted velocity

#include "createPhi1.H"
#include "createPhi2.H"
