// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <boost/filesystem.hpp>
// #include <sbml/SBMLTypes.h>
// #include <sbml/SBMLReader.h>
// SBMLReader reader;

#include <rr/rrRoadRunner.h>

#include <2DTissue.h>

static std::string getSBMLString(){
    return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version2/core\" level=\"3\" version=\"2\">\n"
           "  <model metaid=\"_case00001\" id=\"case00001\" name=\"case00001\" timeUnits=\"time\">\n"
           "    <listOfUnitDefinitions>\n"
           "      <unitDefinition id=\"volume\">\n"
           "        <listOfUnits>\n"
           "          <unit kind=\"litre\" exponent=\"1\" scale=\"0\" multiplier=\"1\"/>\n"
           "        </listOfUnits>\n"
           "      </unitDefinition>\n"
           "      <unitDefinition id=\"substance\">\n"
           "        <listOfUnits>\n"
           "          <unit kind=\"mole\" exponent=\"1\" scale=\"0\" multiplier=\"1\"/>\n"
           "        </listOfUnits>\n"
           "      </unitDefinition>\n"
           "      <unitDefinition id=\"time\">\n"
           "        <listOfUnits>\n"
           "          <unit kind=\"second\" exponent=\"1\" scale=\"0\" multiplier=\"1\"/>\n"
           "        </listOfUnits>\n"
           "      </unitDefinition>\n"
           "    </listOfUnitDefinitions>\n"
           "    <listOfCompartments>\n"
           "      <compartment id=\"compartment\" name=\"compartment\" spatialDimensions=\"3\" size=\"1\" units=\"volume\" constant=\"true\"/>\n"
           "    </listOfCompartments>\n"
           "    <listOfSpecies>\n"
           "      <species id=\"S1\" name=\"S1\" compartment=\"compartment\" initialAmount=\"0.00015\" substanceUnits=\"substance\" hasOnlySubstanceUnits=\"false\" boundaryCondition=\"false\" constant=\"false\"/>\n"
           "      <species id=\"S2\" name=\"S2\" compartment=\"compartment\" initialAmount=\"0\" substanceUnits=\"substance\" hasOnlySubstanceUnits=\"false\" boundaryCondition=\"false\" constant=\"false\"/>\n"
           "    </listOfSpecies>\n"
           "    <listOfParameters>\n"
           "      <parameter id=\"k1\" name=\"k1\" value=\"1\" constant=\"true\"/>\n"
           "    </listOfParameters>\n"
           "    <listOfReactions>\n"
           "      <reaction id=\"reaction1\" name=\"reaction1\" reversible=\"false\">\n"
           "        <listOfReactants>\n"
           "          <speciesReference species=\"S1\" stoichiometry=\"1\" constant=\"true\"/>\n"
           "        </listOfReactants>\n"
           "        <listOfProducts>\n"
           "          <speciesReference species=\"S2\" stoichiometry=\"1\" constant=\"true\"/>\n"
           "        </listOfProducts>\n"
           "        <kineticLaw>\n"
           "          <math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n"
           "            <apply>\n"
           "              <times/>\n"
           "              <ci> compartment </ci>\n"
           "              <ci> k1 </ci>\n"
           "              <ci> S1 </ci>\n"
           "            </apply>\n"
           "          </math>\n"
           "        </kineticLaw>\n"
           "      </reaction>\n"
           "    </listOfReactions>\n"
           "  </model>\n"
           "</sbml>";
}

const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

int main()
{
    int step_count = 30;
    bool save_data = false;
    bool particle_innenleben = false;

    // Path to the 3D mesh file
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";

    rr::RoadRunner roadRunner(getSBMLString());
    std::string filename = PROJECT_PATH.string() + "/sbml-model/BIOMD0000000613_url.xml";
    // SBMLDocument* document = reader.readSBMLFromFile(filename.c_str());

    std::cout << *roadRunner.simulate(0, 1, 11) << std::endl;

    return 0;
}

