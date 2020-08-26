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

Application
    ZMparticleTracks

Description
    Generates a VTK file of particle tracks for cases that were computed using
    a tracked-parcel-type cloud.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "writer.H"
#include "passiveParticleCloud.H"
#include "particleTracksTemplates.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Scanning times to determine track data for cloud " << cloudName
        << nl << endl;

    // max number of originating procs number
    label maxNProcs = 1;
    // the max Id number + 1 for each original processor
    labelList numIds(maxNProcs, 0);
    // allCloud is the combination of different time
    PtrList<passiveParticleCloud> allCloud(0);
    PtrList<IOobjectList> cloudObjsList(0);
	
	label nCloud = 0;
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "    Reading the dictionary of cloudObjs" << endl;

		fileName cloudDir = cloud::prefix/cloudName;
		fileName lagrangianDir = runTime.path()/runTime.timeName()/cloudDir;
        // check the folder in the time step,
        // if true number of cloud nCloud ++;
		if (!isDir(lagrangianDir) || lagrangianDir.empty())
        {
            Info<< "    The dictionary " << lagrangianDir
                << " is empty, jump to next time" << endl;
        }
        else
		{
			cloudObjsList.append
			(
				new IOobjectList
				(
					mesh,
					runTime.timeName(),
					cloud::prefix/cloudName
				)
			);
			allCloud.append
			(
				new passiveParticleCloud
				(
					mesh,
					cloudName
				)
			);
			
			Info<< "    Reading particle positions and origId, origProc" << endl;
			
			passiveParticleCloud& myCloud = allCloud[nCloud];

			Info<< "    Read " << myCloud.size() << " particles" << endl;

			forAllConstIter(passiveParticleCloud, myCloud, iter)
			{
				label origId = iter().origId();
				label origProc = iter().origProc();
				// origProc is starting from 0
				if (origProc >= maxNProcs)
				{
					maxNProcs = origProc+1;
					numIds.setSize(origProc+1, 0);
				}
				numIds[origProc] = max(numIds[origProc], origId+1);
			}
			nCloud ++;
		}
    }

    Info<< nl
		<< "Detected particles originating from " << maxNProcs
        << " processors." << nl << endl;

    Info<< nl << "Particle statistics:" << endl;
    forAll(numIds, proci)
    {
        Info<< "    Found " << numIds[proci] << " particles originating"
            << " from processor " << proci << endl;
    }
    Info<< nl << endl;

    // Calculate starting ids for particles on each processor
    List<label> startIds(maxNProcs, 0);
    for (label proci = 0; proci < maxNProcs-1; proci++)
    {
        startIds[proci+1] += startIds[proci] + numIds[proci];
    }

    // The total number of particles originating from all the processors
    label nParticles = sum(numIds);
    //label nParticles = startIds.last() + numIds[maxNProcs-1];

    // Number of tracks to generate, sampleFrequency is the averaged number of particle
    label nTracks = nParticles/sampleFrequency;
    if (nTracks == 0)
    {
        Info<< "\n    No track data for writting" << endl;
        return 0;
    }

    Info<< "\nGenerating " << nTracks << " particle tracks for cloud "
        << cloudName << nl << endl;

    // Storage the particle index for all sample point in the track lines
    // The first dimension is trackI, second is sampleI
    //List<DynamicList<label>> cloudMap(nTracks);
    //List<DynamicList<label>> particleMap(nTracks);

    // Storage for all particle positions
    List<DynamicList<vector>> trackedPositions(nTracks);
    // Storage the globalParticleI index at
    // all the sample point sampleI in the track lines trackI
    // It is used to write the VTK LINES connections
    List<DynamicList<label>> globalParticleMap(nTracks);

    // A bool list for checking the particle
    // If the particle is tracked
    // sampleParticleMap[iCloud][iParticle] is ture
    List<List<bool>> sampleParticleMap(allCloud.size());

    label globalParticleI = 0;
    forAll(allCloud, iCloud)
    {
        passiveParticleCloud& myCloud = allCloud[iCloud];
        sampleParticleMap[iCloud].setSize(myCloud.size(), false);

        label iParticle = 0;
        forAllConstIter(passiveParticleCloud, myCloud, iter)
        {
            label origId = iter().origId();
            label origProc = iter().origProc();
            label globalId = startIds[origProc] + origId;

            point positions = iter().position();

            if (globalId % sampleFrequency == 0)
            {
                label trackI = globalId/sampleFrequency;
                label trackLength = trackedPositions[trackI].size();
                if (trackLength < maxPositions)
                {
                    // particleMap and cloudMap have a two dimensions
                    // cloudMap[trackI][sampleI] = iCloud;
                    // particleMap[trackI][sampleI] = iParticle;
                    //cloudMap[trackI].append( iCloud );
                    //particleMap[trackI].append( iParticle );

                    trackedPositions[trackI].append( positions );
                    globalParticleMap[trackI].append( globalParticleI );
                    sampleParticleMap[iCloud][iParticle] = true;

                    globalParticleI ++;
                }
            }
            iParticle++;
        }
    }


    //- start to write the data
    if (Pstream::master())
    {
        Info<< "\n    Generating " << nTracks << " tracks for cloud"
            << cloudName << nl << endl;

        fileName vtkPath(runTime.path()/"VTK");
        mkDir(vtkPath);

        OFstream os(vtkPath/"particleTracks.vtk");

        Info<< "\n    Writing particle tracks to " << os.name() << endl;

        os  << "# vtk DataFile Version 2.0" << nl
            << "particleTracks" << nl
            << "ASCII" << nl
            << "DATASET POLYDATA" << nl;

        PtrList<coordSet> tracks(nTracks);
        forAll(tracks, trackI)
        {
            tracks.set
            (
                trackI,
                new coordSet
                (
                    "track" + Foam::name(trackI),
                    "distance"
                )
            );
            // a list of points for the trackI line
            tracks[trackI].transfer
            (
                trackedPositions[trackI]
            );
        }

        label nPoints = 0;
        forAll(tracks, i)
        {
            nPoints += tracks[i].size();
        }

        Info<< "\n    Writing points" << endl;

        os  << "POINTS " << nPoints << " float" << nl;

        forAll(allCloud, iCloud)
        {
            passiveParticleCloud& myCloud = allCloud[iCloud];
            label iParticle = 0;
            forAllConstIter(passiveParticleCloud, myCloud, iter)
            {
                const vector& pt = iter().position();
                if( sampleParticleMap[iCloud][iParticle] )
                {
                    os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
                }
                iParticle++;
            }
        }

        // Write track (line) connectivity to file
        Info<< "\n    Writing track lines" << endl;

        os  << "LINES " << nTracks << ' ' << nPoints+nTracks << nl;

        forAll(tracks, trackI)
        {
            os  << globalParticleMap[trackI].size();
            forAll(globalParticleMap[trackI], sampleI)
            {
                os  << ' ' << globalParticleMap[trackI][sampleI];
            }
            os << nl;
        }

		// Extract the fields name list		
		const List<word> labelFieldNames = validFields<label>
        (
            selectedLagrangianFieldNames,
            cloudObjsList
        );
		const List<word> scalarFieldNames = validFields<scalar>
        (
            selectedLagrangianFieldNames,
            cloudObjsList
        );
		const List<word> vectorFieldNames = validFields<vector>
        (
            selectedLagrangianFieldNames,
            cloudObjsList
        );
		const List<word> sphericalTensorFieldNames = validFields<sphericalTensor>
        (
            selectedLagrangianFieldNames,
            cloudObjsList
        );
		const List<word> symmTensorFieldNames = validFields<symmTensor>
        (
            selectedLagrangianFieldNames,
            cloudObjsList
        );
		const List<word> tensorFieldNames = validFields<tensor>
        (
            selectedLagrangianFieldNames,
            cloudObjsList
        );

		const label nFields = labelFieldNames.size()
							+ scalarFieldNames.size()
							+ vectorFieldNames.size()
							+ sphericalTensorFieldNames.size()
							+ symmTensorFieldNames.size()
							+ tensorFieldNames.size();

        os  << "POINT_DATA " << nPoints << nl
            << "FIELD attributes " << nFields << nl;

        Info<< "\n    Processing fields" << nl << endl;

		Info<< "    Processing label fields" << endl;
		forAll(labelFieldNames, iName)
		{
			const word& fieldName = labelFieldNames[iName];
			os  << nl << fieldName << ' ' << pTraits<label>::nComponents
				<< ' ' << nPoints << " float" << nl;
			processField<label>(os, fieldName, sampleParticleMap, cloudObjsList);
		}

		Info<< "    Processing scalar fields" << endl;
		forAll(scalarFieldNames, iName)
		{
			const word& fieldName = scalarFieldNames[iName];
			os  << nl << fieldName << ' ' << pTraits<scalar>::nComponents
				<< ' ' << nPoints << " float" << nl;
			processField<scalar>(os, fieldName, sampleParticleMap, cloudObjsList);
		}

		Info<< "    Processing vector fields" << endl;
		forAll(vectorFieldNames, iName)
		{
			const word& fieldName = vectorFieldNames[iName];
			os  << nl << fieldName << ' ' << pTraits<vector>::nComponents
				<< ' ' << nPoints << " float" << nl;
			processField<vector>(os, fieldName, sampleParticleMap, cloudObjsList);
		}

		Info<< "    Processing sphericalTensor fields" << endl;
		forAll(sphericalTensorFieldNames, iName)
		{
			const word& fieldName = sphericalTensorFieldNames[iName];
			os  << nl << fieldName << ' ' << pTraits<sphericalTensor>::nComponents
				<< ' ' << nPoints << " float" << nl;
			processField<sphericalTensor>(os, fieldName, sampleParticleMap, cloudObjsList);
		}

		Info<< "    Processing symmTensor fields" << endl;
		forAll(symmTensorFieldNames, iName)
		{
			const word& fieldName = symmTensorFieldNames[iName];
			os  << nl << fieldName << ' ' << pTraits<symmTensor>::nComponents
				<< ' ' << nPoints << " float" << nl;
			processField<symmTensor>(os, fieldName, sampleParticleMap, cloudObjsList);
		}

		Info<< "    Processing tensor fields" << endl;
		forAll(tensorFieldNames, iName)
		{
			const word& fieldName = tensorFieldNames[iName];
			os  << nl << fieldName << ' ' << pTraits<tensor>::nComponents
				<< ' ' << nPoints << " float" << nl;
			processField<tensor>(os, fieldName, sampleParticleMap, cloudObjsList);
		}

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
