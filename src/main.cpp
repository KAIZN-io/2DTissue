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

    // Get the list of species
    ListOfSpecies* speciesList = model->getListOfSpecies();

    for (unsigned int i = 0; i < speciesList->size(); ++i)
    {
        Species* species = dynamic_cast<Species*>(speciesList->get(i));
        std::cout << "Species: " << species->getId() << std::endl;
    }

    // Get the list of reactions
    ListOfReactions* reactionList = model->getListOfReactions();

    for (unsigned int i = 0; i < reactionList->size(); ++i)
    {
        Reaction* reaction = dynamic_cast<Reaction*>(reactionList->get(i));
        std::cout << "Reaction: " << reaction->getId() << std::endl;

        // Get the kinetic law of the reaction
        const KineticLaw* kineticLaw = reaction->getKineticLaw();
        if (kineticLaw != NULL) {
            std::cout << "Kinetic Law: " << kineticLaw->getFormula() << std::endl;
        }
    }

    // Deallocate memory
    delete document;

    return 0;
}

