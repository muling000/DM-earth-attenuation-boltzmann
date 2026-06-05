find_package(Doxygen REQUIRED dot)
find_package(Sphinx REQUIRED)
If(NOT BUILD_DOCS_SKIP_PDF)
  find_package(LATEX COMPONENTS PDFLATEX)
  find_package(Latexmk)
endif()

include(GNUInstallDirs)

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/darkprop)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_HTML_OUTPUT_DIR ${SPHINX_BUILD}/html)
set(SPHINX_INDEX_FILE ${SPHINX_HTML_OUTPUT_DIR}/index.html)
set(SPHINX_CONF
  -Dbreathe_default_project=${PROJECT_NAME}
  -Dbreathe_projects.${PROJECT_NAME}=${DOXYGEN_OUTPUT_DIR}/xml
  -Dproject=${PROJECT_NAME}
  -Dversion=${PROJECT_VERSION}
  -Drelease=${PROJECT_VERSION}
)
file(GLOB_RECURSE DOC_SOURCES ${SPHINX_SOURCE}/*.rst)
file(GLOB_RECURSE DARKPROP_HEADERS ${DARKPROP_HEADER_DIR}/*.hpp)

set(PARTICLE_GRAPH ${SPHINX_SOURCE}/_static/images/inherit_graph_particle.png)
set(PARTICLE_GRAPH_IN ${DOXYGEN_OUTPUT_DIR}/html/inherit_graph_5.png)
set(MEDIUM_GRAPH ${SPHINX_SOURCE}/_static/images/inherit_graph_medium.png)
set(MEDIUM_GRAPH_IN ${DOXYGEN_OUTPUT_DIR}/html/inherit_graph_3.png)
set(INJECTOR_GRAPH ${SPHINX_SOURCE}/_static/images/inherit_graph_injector.png)
set(INJECTOR_GRAPH_IN ${DOXYGEN_OUTPUT_DIR}/html/inherit_graph_1.png)
set(INHERIT_GRAPHS ${PARTICLE_GRAPH} ${MEDIUM_GRAPH} ${INJECTOR_GRAPH})
set(INHERIT_GRAPHS_IN ${PARTICLE_GRAPH_IN} ${MEDIUM_GRAPH_IN} ${INJECTOR_GRAPH_IN})

# Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

add_custom_command(
  OUTPUT ${DOXYGEN_INDEX_FILE} ${INHERIT_GRAPHS_IN}
  DEPENDS ${DARKPROP_HEADERS}
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
  MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
  COMMENT "Generating doxygen docs"
  VERBATIM
)

add_custom_target(doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE} ${INHERIT_GRAPHS_IN})

add_custom_command(
  OUTPUT ${INHERIT_GRAPHS}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PARTICLE_GRAPH_IN} ${PARTICLE_GRAPH}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MEDIUM_GRAPH_IN} ${MEDIUM_GRAPH}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${INJECTOR_GRAPH_IN} ${INJECTOR_GRAPH}
  DEPENDS doxygen
  MAIN_DEPENDENCY ${INHERIT_GRAPHS_IN}
  COMMENT "Copying generated image to Sphinx _static directory"
)

add_custom_command(
  OUTPUT ${SPHINX_INDEX_FILE}
  COMMAND ${SPHINX_EXECUTABLE} -M html ${SPHINX_SOURCE} ${SPHINX_BUILD} ${SPHINX_CONF}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${PROJECT_SOURCE_DIR}/README.md ${DOXYGEN_INDEX_FILE} ${INHERIT_GRAPHS}
  MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
  COMMENT "Generating html documentation with Sphinx"
)

add_custom_target(htmldoc DEPENDS ${SPHINX_INDEX_FILE} VERBATIM)
if(LATEX_PDFLATEX_FOUND AND Latexmk_FOUND)
  set(SPHINX_LATEX_OUTPUT_DIR ${SPHINX_BUILD}/latex)
  set(PDF_MANUAL_NAME darkprop-manual)
  set(SPHINX_PDF_FILE ${SPHINX_LATEX_OUTPUT_DIR}/${PDF_MANUAL_NAME}.pdf)
  add_custom_command(
    OUTPUT ${SPHINX_PDF_FILE}
    COMMAND ${SPHINX_EXECUTABLE} -M latexpdf ${SPHINX_SOURCE} ${SPHINX_BUILD} ${SPHINX_CONF}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${PROJECT_SOURCE_DIR}/README.md ${DOXYGEN_INDEX_FILE} ${INHERIT_GRAPHS}
    MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
    COMMENT "Generating pdf documentation with Sphinx"
  )
  add_custom_target(pdfdoc DEPENDS ${SPHINX_PDF_FILE} VERBATIM)
  install(FILES ${SPHINX_PDF_FILE} DESTINATION ${CMAKE_INSTALL_DOCDIR})
  add_custom_target(docs ALL DEPENDS htmldoc pdfdoc VERBATIM)
else()
  add_custom_target(docs ALL DEPENDS htmldoc VERBATIM)
  If(BUILD_DOCS_SKIP_PDF)
    message(INFO "PDF generation is skipped.")
  else()
    message(WARNING "No PdfLaTeX or Latexmk found. PDF output not available.")
  endif()
endif()

add_custom_command(
  TARGET docs POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E remove ${INHERIT_GRAPHS}
)

install(DIRECTORY ${SPHINX_HTML_OUTPUT_DIR} DESTINATION ${CMAKE_INSTALL_DOCDIR})
