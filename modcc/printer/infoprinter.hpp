#pragma once

#include <string>

#include "module.hpp"
#include "printer/printeropt.hpp"

// Build header file comprising mechanism metadata
// and declarations of backend-specific mechanism implementations.

std::string build_info_header(const Module& m, const printer_options& opt, bool cpu=false, bool gpu=false);
