#include <iostream>
#include <filesystem>
#include <string>

#include <utilities/paths_lib.h>

namespace fs = std::filesystem;


std::string get_full_path(const std::string& project_relative_path_str) {
    fs::path first_path = __FILE__;
    fs::path second_path{project_relative_path_str};
   
    std::string project_folder_name = "2DTissue";
    fs::path full_path;

    // Iterate through the first path's components until the project folder is found
    for (const auto& component : first_path) {
        if (component == project_folder_name) {
            break;
        }
        full_path /= component;
    }

    // Append the second path to the new path
    full_path /= second_path;

    // Normalize the path to remove any redundant or repeated components
    full_path = fs::canonical(full_path);

    // Return the path as a string
    return full_path.string();
}

std::string get_repo_root() {
    fs::path first_path = __FILE__;
    std::string project_folder_name = "2DTissue";
    std::filesystem::path repo_root;

    // Iterate through the full path's components until the project folder is found
    for (const auto& component : first_path) {
        if (component == project_folder_name) {
            break;
        }
        repo_root /= component;
    }

    // Return the path as a string
    return repo_root.string();
}