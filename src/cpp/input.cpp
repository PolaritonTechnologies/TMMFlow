#include "input.h"
#include <set>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include "spline.h"

using json = nlohmann::json;

CalculationInfo loadCalculationInfo(const std::string &json_text)
{

    json calculation_order = json::parse(json_text);
    std::vector<std::string> structure_materials = calculation_order["structure_materials"];
    std::vector<bool> incoherent = calculation_order["incoherent"];
    std::vector<double> structure_thicknesses = calculation_order["structure_thicknesses"];

    CalculationInfo calculation_info(
        calculation_order["calculation_type"],
        calculation_order["polarAngleMin"],
        calculation_order["polarAngleMax"],
        calculation_order["polarAngleStep"],
        calculation_order["azimAngleMin"],
        calculation_order["azimAngleMax"],
        calculation_order["azimAngleStep"],
        calculation_order["wavelengthMin"],
        calculation_order["wavelengthMax"],
        calculation_order["wavelengthStep"],
        calculation_order["polarization"],
        structure_materials,
        incoherent,
        structure_thicknesses,
        calculation_order["incident_medium"],
        calculation_order["exit_medium"]);

    return calculation_info;
}

std::vector<std::string> getUniqueMembers(const std::vector<std::string> &inputArray)
{
    std::set<std::string> uniqueSet(inputArray.begin(), inputArray.end());
    std::vector<std::string> uniqueVector(uniqueSet.begin(), uniqueSet.end());
    return uniqueVector;
}

std::filesystem::path getBasePath()
{
    auto currentPath = std::filesystem::current_path();
    // Check if current execution is from /TMM/src
    if (currentPath.filename() == "src")
    {
        return currentPath.parent_path(); // Move up to /TMM
    }
    // Otherwise, assume execution from /TMM or handle other cases as needed
    return currentPath;
}

std::filesystem::path getDatabasePath(const std::string &filename)
{
    auto basePath = getBasePath();
    return basePath / "instance" / filename;
}
