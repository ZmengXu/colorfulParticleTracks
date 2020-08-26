IOdictionary propsDict
(
    IOobject
    (
        "particleTrackProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED
    )
);

const word cloudName(propsDict.lookup("cloudName"));

List<word> selectedLagrangianFieldNames(propsDict.lookup("fields"));

label sampleFrequency(readLabel(propsDict.lookup("sampleFrequency")));

label maxPositions(readLabel(propsDict.lookup("maxPositions")));
