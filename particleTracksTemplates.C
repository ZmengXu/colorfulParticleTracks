/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "particleTracksTemplates.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::fieldOk(const IOobjectList& cloudObjs, const word& name)
{
    IOobjectList objects(cloudObjs.lookupClass(IOField<Type>::typeName));

    return (objects.lookup(name) != nullptr);
}


template<class Type>
Foam::List<Foam::word>
Foam::validFields
(
    const List<word>& fieldNames,
    const PtrList<IOobjectList>& cloudObjsList
)
{
    List<word> validFieldNames(0);
    
    forAll(fieldNames, iName)
    {
        const word& fieldName = fieldNames[iName];
        bool fieldNameOK(true);

        forAll(cloudObjsList, iCloud)
        {
            const IOobjectList& cloudObjs = cloudObjsList[iCloud];
            fieldNameOK = fieldNameOK && fieldOk<Type>(cloudObjs, fieldName);
        }
        
        if (fieldNameOK)
        {
            validFieldNames.append(fieldName);
        }
    }
    return validFieldNames;
}


template<class Type>
void Foam::readField
(
    List<Type>& values,
    const word& fieldName,
    const IOobjectList& cloudObjs
)
{
    IOobjectList objects(cloudObjs.lookupClass(IOField<Type>::typeName));

    const IOobject* obj = objects.lookup(fieldName);
    if (obj != nullptr)
    {
        IOField<Type> newField(*obj);
        values = newField;
    }
    else
    {
        FatalErrorInFunction
            << "Unable to read field " << fieldName
            << abort(FatalError);
    }
}


template<class Type>
void Foam::writeVTK(OFstream& os, const label& value)
{
    os  << value;
}


template<class Type>
void Foam::writeVTK(OFstream& os, const scalar& value)
{
    os  << value;
}


template<class Type>
void Foam::writeVTK(OFstream& os, const Type& value)
{
    os  << value.component(0);
    for (label i=1; i<pTraits<Type>::nComponents; i++)
    {
        os  << ' ' << value.component(i);
    }
}


template<class Type>
void Foam::processField
(
    OFstream& os,
    const word& validFieldName,
    const List<List<bool>>& sampleParticleMap,
    const PtrList<IOobjectList>& cloudObjsList
)
{
    const label step = max(floor(8/pTraits<Type>::nComponents), 1);

	forAll(cloudObjsList, iCloud)
	{
		const IOobjectList& cloudObjs = cloudObjsList[iCloud];
		const IOobject* obj = cloudObjs.lookup(validFieldName);

        Info<< "        reading field " << validFieldName
			<< " in the time path of " << obj->instance()
			<< endl;
			
		List<Type> values;
		readField<Type>
		(
			values,
			validFieldName,
			cloudObjs
		);

		Info<< "        writing field " << validFieldName << endl;

		label nData = values.size() - 1;
		forAll( values, iParticle )
		{
			if( sampleParticleMap[iCloud][iParticle] )
			{
				writeVTK<Type>(os, values[iParticle]);
				if
				(
					((iParticle + 1) % step == 0)
					||
					(iParticle == nData)
				)
				{
					os  << nl;
				}
				else
				{
					os  << ' ';
				}
			}
		}
	}
}

// ************************************************************************* //
