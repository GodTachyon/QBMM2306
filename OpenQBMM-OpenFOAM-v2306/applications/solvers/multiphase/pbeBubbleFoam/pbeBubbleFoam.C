/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Copyright (C) 2015-2020 Alberto Passalacqua
-------------------------------------------------------------------------------
2016-06-30 Alberto Passalacqua: Implemented quadrature-based bubbly flow solver.
                                Implementation performed by J. C. Heylmun
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    pbeBubbleFoam

Description
    Solver for polydisperse gas-liquid flows using a quadrature-based moment
    method to describe the gas phase.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "twoPhaseSystemPbeBubble.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for polydisperse gas-liquid flows using a quadrature-based "
        "moment method to describe the gas phase."
    );

    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    Switch implicitPhasePressure
    (
        mesh.solverDict(alpha1.name()).lookupOrDefault<Switch>
        (
            "implicitPhasePressure", false
        )
    );

    #include "pU/createDDtU.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

      {

            // Solve for mean phase velocities and gas volume fraction
            while (pimple.loop())
            {
                fluid.solve();
                fluid.correct();

                #include "contErrs.H"
                #include "pU/DDtU.H"
                #include "pU/UEqns.H"

                #include "pU/pEqn.H"

                if (pimple.turbCorr())
                {
                    // Transport moments with mean gas velocity
                    fluid.averageTransport();
                    fluid.correctTurbulence();
                }
            }
        }
        /*
        // [NEW ADDITION] Print moments and dimensions for diagnostics
        {
            const polydispersePhaseModel& pbePhase1 =
                refCast<const polydispersePhaseModel>(phase1);

            const volScalarMomentFieldSet& moments1 =
                pbePhase1.quadrature().moments();

            Info<< "Phase1 moments:" << endl;
            forAll(moments1, mi)
            {
                const volScalarField& m = moments1[mi];
                Info<< "    " << m.name()
                    << "  dims: " << m.dimensions()
                    << "  min: " << Foam::min(m).value()
                    << "  max: " << Foam::max(m).value()
                    << "  mean: "
                    << m.weightedAverage(mesh.V()).value()
                    << endl;
            }
        }
        */
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
