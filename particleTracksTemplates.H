/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#ifndef particleTracksTemplates_H
#define particleTracksTemplates_H

#include "OFstream.H"
#include "SortableList.H"
#include "IOobjectList.H"
#include "PtrList.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool fieldOk(const IOobjectList& cloudObjs, const word& name);

template<class Type>
void writeVTK(OFstream& os, const label& value);

template<class Type>
void writeVTK(OFstream& os, const scalar& value);

template<class Type>
void writeVTK(OFstream& os, const Type& value);

template<class Type>
List<word> validFields
(
    const List<word>& fieldNames,
    const PtrList<IOobjectList>& cloudObjsList
);

template<class Type>
void readField
(
    List<Type>& values,
    const word& fieldName,
    const IOobjectList& cloudObjs
);

template<class Type>
void processField
(
    OFstream& os,
    const word& fieldName,
    const List<List<bool>>& sampleParticleMap,
    const PtrList<IOobjectList>& cloudObjsList
);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "particleTracksTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
