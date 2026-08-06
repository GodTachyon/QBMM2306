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
    const auto& quadrature = pbe_->quadrature();

    d_ = dimensionedScalar("zero", dimLength, 0.0);
    // 3rd diameter moment. Have to check if this is correct in dimensions
    volScalarField d32Num
    (
        IOobject
        (
            "d32Num",
            fluid_.mesh().time().timeName(),
            fluid_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fluid_.mesh(),
        dimensionedScalar
        (
            "zero",
            quadrature.nodes()[0].primaryWeight().dimensions()*pow3(dimLength),
            0.0
        )
    );
    // 2nd diameter moment. Have to check dimensions.
    volScalarField d32Den
    (
        IOobject
        (
            "d32Den",
            fluid_.mesh().time().timeName(),
            fluid_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fluid_.mesh(),
        dimensionedScalar
        (
            "zero",
            quadrature.nodes()[0].primaryWeight().dimensions()*sqr(dimLength),
            0.0
        )
    );

    forAll(quadrature.nodes(), nodei)
    {
        const auto& node = quadrature.nodes()[nodei];

        alphas_[nodei] =
            node.primaryWeight()*node.primaryAbscissae()[0]/rho();
        alphas_[nodei].max(0);
        alphas_[nodei].min(1);

        ds_[nodei] =
            Foam::min
            (
                Foam::max
                (
                    Foam::pow
                    (
                        node.primaryAbscissae()[0]*6.0
                       /(rho()*Foam::constant::mathematical::pi)
                      + dimensionedScalar("smallVolume", dimVolume, SMALL),
                        1.0/3.0
                    ),
                    minD_
                ),
                maxD_
            );

        d_ += alphas_[nodei]*ds_[nodei];
        d32Num += node.primaryWeight()*pow3(ds_[nodei]);
        d32Den += node.primaryWeight()*sqr(ds_[nodei]);
    }

    d_.max(minD_);

    d32_ =
        d32Num
       /Foam::max(d32Den, dimensionedScalar("dSmall", d32Den.dimensions(), SMALL));
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
