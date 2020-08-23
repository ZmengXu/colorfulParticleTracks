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
#include "passiveParticleCloud.H"
#include "writer.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::addOption
    (
        "lagrangianFields",
        "list",
        "specify a list of lagrangian fields to be reconstructed. Eg, '(d T)' -"
        "only scalar field is supported currently, "
        "positions always included."
    );
    #include "setRootCase.H"

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fileName vtkPath(runTime.path()/"VTK");
    mkDir(vtkPath);

    Info<< "Scanning times to determine track data for cloud " << cloudName
        << nl << endl;

    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "    Reading particle positions" << endl;
        passiveParticleCloud myCloud(mesh, cloudName);

        Info<< "    Read " << returnReduce(myCloud.size(), sumOp<label>())
            << " particles" << endl;

        forAllConstIter(passiveParticleCloud, myCloud, iter)
        {
            label origId = iter().origId();
            label origProc = iter().origProc();

            if (origProc >= maxIds.size())
            {
                maxIds.setSize(origProc+1, -1);
            }

            maxIds[origProc] = max(maxIds[origProc], origId);
        }
    }

    label maxNProcs = returnReduce(maxIds.size(), maxOp<label>());

    Info<< "Detected particles originating from " << maxNProcs
        << " processors." << nl << endl;

    maxIds.setSize(maxNProcs, -1);

    Pstream::listCombineGather(maxIds, maxEqOp<label>());
    Pstream::listCombineScatter(maxIds);

    labelList numIds = maxIds + 1;

    Info<< nl << "Particle statistics:" << endl;
    forAll(maxIds, proci)
    {
        Info<< "    Found " << numIds[proci] << " particles originating"
            << " from processor " << proci << endl;
    }
    Info<< nl << endl;


    // Calculate starting ids for particles on each processor
    List<label> startIds(numIds.size(), 0);
    for (label i = 0; i < numIds.size()-1; i++)
    {
        startIds[i+1] += startIds[i] + numIds[i];
    }
    label nParticle = startIds.last() + numIds[startIds.size()-1];


    // Number of tracks to generate
    label nTracks = nParticle/sampleFrequency;

    // Storage for all particle tracks
    List<DynamicList<vector>> allTracks(nTracks);


    List<word> selectedLagrangianFields;
    if (args.optionFound("lagrangianFields"))
    {
        args.optionLookup("lagrangianFields")() >> selectedLagrangianFields;
    }

    Info<< "\nLooking for the lagrangianFields " 
        << selectedLagrangianFields << nl << endl;

    // Storage for all fields tracks
    List<List<DynamicList<scalar>>> allTrackFields(selectedLagrangianFields.size());
    forAll(allTrackFields, fieldI)
    {
        allTrackFields[fieldI].setSize
        (
            nTracks
        );
    }

    Info<< "\nGenerating " << nTracks << " particle tracks for cloud "
        << cloudName << nl << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        List<pointField> allPositions(Pstream::nProcs());
        List<labelField> allOrigIds(Pstream::nProcs());
        List<labelField> allOrigProcs(Pstream::nProcs());

        // Read particles. Will be size 0 if no particles.
        Info<< "    Reading particle positions" << endl;
        passiveParticleCloud myCloud(mesh, cloudName);

        // Collect the track data on all processors that have positions
        allPositions[Pstream::myProcNo()].setSize
        (
            myCloud.size(),
            point::zero
        );
        allOrigIds[Pstream::myProcNo()].setSize(myCloud.size(), 0);
        allOrigProcs[Pstream::myProcNo()].setSize(myCloud.size(), 0);

        label i = 0;
        forAllConstIter(passiveParticleCloud, myCloud, iter)
        {
            allPositions[Pstream::myProcNo()][i] = iter().position();
            allOrigIds[Pstream::myProcNo()][i] = iter().origId();
            allOrigProcs[Pstream::myProcNo()][i] = iter().origProc();
            i++;
        }

        // Collect the track data on the master processor
        Pstream::gatherList(allPositions);
        Pstream::gatherList(allOrigIds);
        Pstream::gatherList(allOrigProcs);
        
        PtrList<List<scalarField>> allLagrangianfields(selectedLagrangianFields.size());
        forAll(selectedLagrangianFields, fieldI)
        {
            allLagrangianfields.set
            (
                fieldI,
                new List<scalarField>
                (
                    Pstream::nProcs()
                )
            );

            const word fieldName = selectedLagrangianFields[fieldI];
            // Check object on local mesh
            IOobject fieldIOobject
            (
                fieldName,
                runTime.timeName(),
                cloud::prefix/cloudName,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (fieldIOobject.typeHeaderOk<IOField<scalar>>(true))
            {
                IOField<scalar> scalarIOField(fieldIOobject);
                Field<scalar>& scalarField(scalarIOField);
                allLagrangianfields[fieldI][Pstream::myProcNo()]=scalarField;
            }

            Pstream::gatherList(allLagrangianfields[fieldI]);
        }


        Info<< "    Constructing tracks" << nl << endl;
        if (Pstream::master())
        {
            forAll(allPositions, proci)
            {
                forAll(allPositions[proci], i)
                {
                    label globalId =
                        startIds[allOrigProcs[proci][i]]
                      + allOrigIds[proci][i];

                    if (globalId % sampleFrequency == 0)
                    {
                        label trackId = globalId/sampleFrequency;
                        if (allTracks[trackId].size() < maxPositions)
                        {
                            allTracks[trackId].append
                            (
                                allPositions[proci][i]
                            );

                            forAll(allLagrangianfields, fieldI)
                            {
                                allTrackFields[fieldI][trackId].append
                                (
                                    allLagrangianfields[fieldI][proci][i]
                                );
                            }
                        }
                    }
                }
            }
        }
    }


    if (Pstream::master())
    {
        PtrList<coordSet> tracks(allTracks.size());
        forAll(allTracks, trackI)
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
            tracks[trackI].transfer(allTracks[trackI]);
        }

        List<List<scalarField>> trackFields(allTrackFields.size());
        forAll(trackFields, fieldI)
        {
            trackFields[fieldI].setSize
            (
                nTracks
            );
            forAll(trackFields[fieldI], trackI)
            {
                trackFields[fieldI][trackI] = allTrackFields[fieldI][trackI];
            }
        }

        autoPtr<writer<scalar>> scalarFormatterPtr = writer<scalar>::New
        (
            setFormat
        );

        // OFstream vtkTracks(vtkPath/"particleTracks.vtk");
        fileName vtkFile
        (
            scalarFormatterPtr().getFileName
            (
                tracks[0],
                wordList(0)
            )
        );

        OFstream vtkTracks
        (
            vtkPath/("particleTracks." + vtkFile.ext())
        );

        Info<< "\nWriting particle tracks in " << setFormat
            << " format to " << vtkTracks.name()
            << nl << endl;

        scalarFormatterPtr().write
        (
            true,
            tracks,
            selectedLagrangianFields,//wordList(trackFields.size()),
            trackFields,//List<List<scalarField>>(0),
            vtkTracks
        );
    }

    return 0;
}


// ************************************************************************* //
