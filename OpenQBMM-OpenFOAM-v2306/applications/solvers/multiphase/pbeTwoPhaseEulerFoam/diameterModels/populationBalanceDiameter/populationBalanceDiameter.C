#include "populationBalanceDiameter.H"
#include "addToRunTimeSelectionTable.H"
#include "populationBalanceModel.H"
#include "univariatePDFTransportModel.H"

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(populationBalance, 0);
    addToRunTimeSelectionTable(diameterModel, populationBalance, dictionary);
}
}

Foam::diameterModels::populationBalance::populationBalance
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    minD_("minD", dimLength, diameterProperties_),
    maxD_("maxD", dimLength, diameterProperties_)
{}

Foam::diameterModels::populationBalance::~populationBalance()
{}

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalance::d() const
{
    const fvMesh& mesh = phase_.U().mesh();

    const HashTable<const Foam::populationBalanceModel*> pbes
    (
        mesh.lookupClass<Foam::populationBalanceModel>()
    );

    if (pbes.empty())
    {
        FatalErrorInFunction
            << "No populationBalanceModel is registered on the mesh."
            << exit(FatalError);
    }

    const auto* pbePtr =
        dynamic_cast<const PDFTransportModels::univariatePDFTransportModel*>
        (
            *pbes.begin()
        );

    if (!pbePtr)
    {
        FatalErrorInFunction
            << "Registered populationBalanceModel is not a "
            << "univariatePDFTransportModel; cannot compute d32."
            << exit(FatalError);
    }

    const volScalarField& m2 = pbePtr->quadrature().moments()[2];
    const volScalarField& m3 = pbePtr->quadrature().moments()[3];

    tmp<volScalarField> td32
    (
        new volScalarField
        (
            IOobject
            (
                "d",
                phase_.U().time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            m3/max(m2, dimensionedScalar("dSmall", m2.dimensions(), SMALL))
        )
    );

    td32.ref().min(maxD_);
    td32.ref().max(minD_);

    return td32;
}

bool Foam::diameterModels::populationBalance::read
(
    const dictionary& diameterProperties
)
{
    diameterModel::read(diameterProperties);

    diameterProperties_.readEntry("minD", minD_);
    diameterProperties_.readEntry("maxD", maxD_);

    return true;
}
