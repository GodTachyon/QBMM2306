/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2022 Alberto Passalacqua
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
    
Class
    Foam::pbeBubblePhaseModel

Description
    Polydisperse phase model that delegates size-moment transport,
    coalescence, and breakup to a stock OpenQBMM
    PDFTransportModels::populationBalanceModels::univariatePopulationBalance
    object, instead of hand-rolling quadrature/kernel/ODE logic locally
    (as polydispersePhaseModel does).

    Built as a second, parallel phaseModel to compare against
    polydispersePhaseModel. This model:
      - Uses scalarQuadratureApproximation (no per-node velocity moments)
        instead of monoKineticQuadratureApproximation, matching the
        single-shared-gas-velocity design.
      - Coalescence/breakup/growth/nucleation come from the stock
        aggregationKernel/breakupKernel/growthModel/nucleationModel
        selection tables via univariatePopulationBalance, rather than
        the phase's own kernel members.

    ASSUMPTIONS requiring verification (see chat discussion):
      - scalarQuadratureApproximation's node/moment accessors mirror
        monoKineticQuadratureApproximation's (primaryWeight(),
        primaryAbscissae()[0], moments()[i]). Not confirmed against
        quadratureApproximations.H / momentFieldSets.H.
      - No corr_/nFluxCorrectors-style consistency enforcement between
        alpha1 and the PBE's own moment[1] is implemented here (unlike
        polydispersePhaseModel). Left out deliberately for this first
        comparison pass.

SourceFiles
    pbeBubblePhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef pbeBubblePhaseModel_H
#define pbeBubblePhaseModel_H

#include "phaseModel.H"
#include "univariatePopulationBalance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class pbeBubblePhaseModel Declaration
\*---------------------------------------------------------------------------*/

class pbeBubblePhaseModel
:
    public phaseModel
{
    // Private data

        //- Population balance properties dictionary
        //  NOTE: schema differs from polydispersePhaseModel's
        //  populationBalanceProperties - this one follows the stock
        //  univariatePopulationBalance format (aggregation/breakup/growth/
        //  nucleation switches, aggregationKernel/breakupKernel/growthModel/
        //  diffusionModel/nucleationModel subdicts, "univariateCoeffs").
        IOdictionary pbeDict_;

        //- Stock univariate population balance model
        autoPtr
        <
            PDFTransportModels::populationBalanceModels
            ::univariatePopulationBalance
        > pbe_;

        //- Number of quadrature nodes
        const label nNodes_;

        //- List of alphas per node
        PtrList<volScalarField> alphas_;

        //- List of diameters per node
        PtrList<volScalarField> ds_;

        //- Maximum diameter
        dimensionedScalar maxD_;

        //- Minimum diameter
        dimensionedScalar minD_;

        //- Sauter mean diameter
        volScalarField d32_;


public:

    //- Runtime type information
    TypeName("pbeBubblePhaseModel");


    // Constructors

        pbeBubblePhaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& phaseProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~pbeBubblePhaseModel();


    // Member Functions

        //- Return the number of nodes
        virtual label nNodes() const
        {
            return nNodes_;
        }

        //- Return alpha field for nodei
        virtual const volScalarField& alphas(const label nodei) const
        {
            return alphas_[nodei];
        }

        //- Return non-constant access to alpha field for nodei
        virtual volScalarField& alphas(const label nodei)
        {
            return alphas_[nodei];
        }

        //- Return the diameter for nodei
        virtual const volScalarField& ds(const label nodei) const
        {
            return ds_[nodei];
        }

        //- Correct the phase properties (alphas_, ds_, d32_) from the
        //  PBE's current quadrature
        virtual void correct();

        //- No-op: single shared gas velocity means there is no relative
        //  flux to transport moments with
        virtual void relativeTransport()
        {
            return;
        }

        //- Advance the population balance with phiGas, then correct()
        virtual void averageTransport(const surfaceScalarField& phiGas);

        //- Read phase properties dictionary
        virtual bool read(const bool readOK);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
