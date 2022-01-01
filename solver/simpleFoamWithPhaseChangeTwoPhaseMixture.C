/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    interPhaseChangeFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids with phase-change
    (e.g. cavitation).  Uses a VOF (volume of fluid) phase-fraction based
    interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    The set of phase-change models provided are designed to simulate cavitation
    but other mechanisms of phase-change are supported within this solver
    framework.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"
#include "myPhaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "newForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", dimMass/dimTime, 0)
            );

            mixture->correct();

            #include "alphaEqn.H"

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();
  //    myPoints.regIOobject::write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    
    /**************************************** commented section  *******************************************/
    /*
     
    // if you want to calculate forces and momentum on a specific patch
    
    wordReList patch_(1); 
	patch_[0]="name of the patch";  // given by user
	const point centreOfRotat_(0,0,0);  // the center around which the torque is calculated (given by user)
	
	dictionary forcesDict;  // create a dictionary to give required data (below)

	forcesDict.add("type", "forces");
	forcesDict.add("patches", patch_);
	forcesDict.add("rhoInf", 1.0);
	forcesDict.add("rho", "rhoInf");
	forcesDict.add("CofR", centreOfRotat_);

	functionObjects::newForces f("newForces", mesh , forcesDict); // create an object from newForces class
																  // to calculate forces and momentum 

	f.calcForcesMoment();   // forces and momentum are calculated
	
	// create a file to write forces and momentum on a specific patch (below) 

	OFstream os("name of the file.txt");   // address of the output file (given by user)
	
	os<<"Fx = "<< f.forceEff().x()<< endl;
	os<<"Fy = "<< f.forceEff().y()<< endl;
	os<<"Fz = "<< f.forceEff().z()<< endl;
	os<<"Mx = "<< f.momentEff().x()<< endl;
	os<<"My = "<< f.momentEff().y()<< endl;
	os<<"Mz = "<< f.momentEff().z()<< endl;     
	             
	*/
	
	/***************************************************************************************/

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
