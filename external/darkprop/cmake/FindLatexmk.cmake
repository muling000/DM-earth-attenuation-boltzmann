find_program(LATEXMK_EXECUTABLE
             NAMES latexmk
             DOC "Path to latexmk executable")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Latexmk
                                  "Failed to find latexmk executable"
                                  LATEXMK_EXECUTABLE)
