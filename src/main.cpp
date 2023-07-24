// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <boost/filesystem.hpp>
#include <sbml/SBMLTypes.h>
#include <sbml/SBMLReader.h>

SBMLReader reader;

#include <2DTissue.h>

const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

int main()
{
    int step_count = 30;
    bool save_data = false;
    bool particle_innenleben = false;

    // Path to the 3D mesh file
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";

    std::string filename = PROJECT_PATH.string() + "/sbml-model/BIOMD0000000613_url.xml";
    SBMLDocument* document = reader.readSBMLFromFile(filename.c_str());

    // If there are errors in the document, print them and stop.
    if (document->getNumErrors() > 0) {
        document->printErrors();
        return 1;
    }

    // Get the model from the document.
    Model* model = document->getModel();

    // Do something with the model, for example print its name
    if(model != NULL) {
        std::cout << "Loaded model: " << model->getId() << std::endl;
    } else {
        std::cout << "No model present." << std::endl;
    }

    // Deallocate memory
    delete document;

    return 0;
}
