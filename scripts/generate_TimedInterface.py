# Generate definition of the TimedInterface C++ class

from pathlib import Path


HEADER = """#include "TimedInterface.hpp"

#include <chrono>

using namespace std;
using namespace discamb;
using namespace chrono;

namespace py = pybind11;

std::vector<Runtime> TimedInterface::get_runtimes() const {
    return mRuntimes;
}
"""

RETURN_TEMPLATE = """
{return_type} TimedInterface::{function_name}({args}) {{
    auto runtime_start = high_resolution_clock::now();
    {return_type} res = PythonInterface::{function_name}({args_values});
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({{std::string("{function_name}"), delta.count()}});
    return res;
}}
"""
VOID_TEMPLATE = """
void TimedInterface::{function_name}({args}) {{
    auto runtime_start = high_resolution_clock::now();
    PythonInterface::{function_name}({args_values});
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({{std::string("{function_name}"), delta.count()}});
}}
"""


def declaration_to_implementation(dec: str) -> str:
    return_type = dec.split(" ")[0].strip()
    function_name = dec[len(return_type) + 1 : dec.index("(")].strip()
    args = dec[dec.index("(") + 1 : dec.index(")")].strip()
    args_values = ", ".join(a.split(" ")[-1] if a else a for a in args.split(", "))

    template = VOID_TEMPLATE if return_type == "void" else RETURN_TEMPLATE
    return template.format(
        return_type=return_type,
        function_name=function_name,
        args=args,
        args_values=args_values,
    )


def main():
    hpp = Path(__file__).parent.parent / "include" / "TimedInterface.hpp"
    cpp = Path(__file__).parent.parent / "src" / "TimedInterface.cpp"

    declarations = []
    found_start = False
    for line in hpp.read_text().split("\n"):
        if "// clang-format" in line:
            found_start ^= True
            continue
        if found_start:
            declarations.append(line.strip())

    with cpp.open("w") as f:
        f.write(HEADER)
        f.write("\n// clang-format off")
        for dec in declarations:
            f.write(declaration_to_implementation(dec))
        f.write("// clang-format on\n")


if __name__ == "__main__":
    main()
