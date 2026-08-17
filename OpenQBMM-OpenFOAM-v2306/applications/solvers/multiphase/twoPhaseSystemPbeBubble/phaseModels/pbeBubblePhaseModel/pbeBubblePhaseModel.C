/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "pbeBubblePhaseModel.H"
#include "twoPhaseSystemPbeBubble.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pbeBubblePhaseModel, 0);
//  addToRunTimeSelectionTable(pbeBubblePhaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pbeBubblePhaseModel::pbeBubblePhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid, phaseProperties, phaseName),
    pbeDict_
    (
        IOobject
        (
            "populationBalanceProperties",
            fluid.mesh().time().constant(),
            fluid.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    pbe_
    (
        new PDFTransportModels::populationBalanceModels
            ::univariatePopulationBalance
            (
                phaseName,
                pbeDict_.subDict("univariateCoeffs"),
                this->phi()
            )
    ),
    nNodes_(pbe_->quadrature().nodes().size()),
    alphas_(nNodes_),
    ds_(nNodes_),
    maxD_("maxD", dimLength, phaseDict_),
    minD_("minD", dimLength, phaseDict_),
    d32_
    (
        IOobject
        (
            "d32",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("d32", dimLength, minD_.value())
    )
{
    this->d_.writeOpt() = IOobject::AUTO_WRITE;

    forAll(alphas_, nodei)
    {
        alphas_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "alpha",
                        IOobject::groupName(phaseModel::name_, Foam::name(nodei))
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("alpha", dimless, 0.0)
            )
        );

        ds_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "d",
                        IOobject::groupName(phaseModel::name_, Foam::name(nodei))
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh(),
                minD_
            )
        );
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pbeBubblePhaseModel::~pbeBubblePhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pbeBubblePhaseModel::correct()
{
    // ASSUMPTION: scalarQuadratureApproximation's node/moment accessors
    // mirror monoKineticQuadratureApproximation's, as used in
    // polydispersePhaseModel::correct(). Unverified - check
    // quadratureApproximations.H / momentFieldSets.H if this fails to
    // compile.
    //
    // This quadrature's abscissa is diameter, not mass (unlike
    // polydispersePhaseModel's monoKineticQuadratureApproximation), so:
    //   - ds_[nodei] is the (bounded) abscissa directly - no mass/rho
    //     conversion.
    //   - Volume fraction comes from moment 3 (alpha = (pi/6)*moment[3]
    //     for spherical bubbles), not moment 1, and needs no /rho() since
    //     d^3 already encodes volume geometrically.
    //   - d32 = moment[3]/moment[2] directly, no per-node reconstruction
    //     needed.
    // ASSUMPTION: this requires the case to track at least moments 0-3.
    const auto& quadrature = pbe_->quadrature();
 
    const scalar volCoeff = Foam::constant::mathematical::pi/6.0;
 
    // Raw (unscaled) volume-fraction estimate from moment 3, used to
    // reconcile alphas_ against the independently-transported alpha1
    // field, mirroring polydispersePhaseModel's scale factor.
    volScalarField scale
    (
        (*this)/Foam::max(volCoeff*quadrature.moments()[3], residualAlpha_)
    );
 
    d_ = dimensionedScalar("zero", dimLength, 0.0);
 
    forAll(quadrature.nodes(), nodei)
    {
        const auto& node = quadrature.nodes()[nodei];
 
        ds_[nodei] =
            Foam::min(Foam::max(node.primaryAbscissae()[0], minD_), maxD_);
 
        alphas_[nodei] =
            volCoeff*node.primaryWeight()*pow3(ds_[nodei])*scale;
        alphas_[nodei].max(0);
        alphas_[nodei].min(1);
 
        d_ += alphas_[nodei]*ds_[nodei];
    }
 
    d_.max(minD_);
 
    d32_ =
        quadrature.moments()[3]
       /Foam::max
        (
            quadrature.moments()[2],
            dimensionedScalar
            (
                "dSmall",
                quadrature.moments()[2].dimensions(),
                SMALL
            )
        );
    d32_ = Foam::min(Foam::max(d32_, minD_), maxD_);
 
    Info<< "d32 (Sauter mean diameter, pbeBubblePhaseModel): min = "
        << Foam::min(d32_).value()
        << " max = " << Foam::max(d32_).value()
        << " mean = " << d32_.weightedAverage(d32_.mesh().V()).value()
        << " m" << endl;
}


void Foam::pbeBubblePhaseModel::averageTransport(const surfaceScalarField& phiGas)
{
    // pbe_ was constructed against this->phi(). twoPhaseSystem::
    // averageTransport() calls phase1_->averageTransport(phase1_->phi()),
    // so phiGas should be the exact same field object. If that call site
    // ever changes, this catches the mismatch instead of silently
    // advecting moments with the wrong flux.
    if (&phiGas != &(this->phi()))
    {
        FatalErrorInFunction
            << "pbeBubblePhaseModel's population balance was constructed "
            << "against a different flux than phiGas passed in here."
            << exit(FatalError);
    }

    Info<< "Transporting moments with univariatePopulationBalance" << endl;
    pbe_->solve();

    correct();
}


bool Foam::pbeBubblePhaseModel::read(const bool readOK)
{
    bool read = false;
    if (readOK)
    {
        maxD_.readIfPresent(phaseDict_);
        minD_.readIfPresent(phaseDict_);
        read = true;
    }

    if (pbe_->readIfModified())
    {
        read = true;
    }

    return read || readOK;
}


// ************************************************************************* //
