
#include <fstream>
#include <iostream>

namespace io
{
    std::ofstream prepare_log_file(const std::filesystem::path& data_dir,
                               std::string& log_path);
} // namespace io

