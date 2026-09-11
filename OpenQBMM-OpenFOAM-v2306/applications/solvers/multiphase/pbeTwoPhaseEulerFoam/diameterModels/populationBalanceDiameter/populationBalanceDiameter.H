/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::diameterModels::constant

Description
    Constant dispersed-phase particle diameter model.

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/
#ifndef diameterModels_populationBalance_H
#define diameterModels_populationBalance_H

#include "diameterModel.H"

namespace Foam
{
namespace diameterModels
{

class populationBalance
:
    public diameterModel
{
    // Private data

        //- Minimum diameter
        dimensionedScalar minD_;

        //- Maximum diameter
        dimensionedScalar maxD_;

public:

    TypeName("populationBalance");

    populationBalance
    (
        const dictionary& diameterProperties,
        const phaseModel& phase
    );

    virtual ~populationBalance();

    virtual tmp<volScalarField> d() const;

    virtual bool read(const dictionary& diameterProperties);
};

} // End namespace diameterModels
} // End namespace Foam

#endif

// ************************************************************************* //
