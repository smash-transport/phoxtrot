# project name
project(phoxtrot)

# minimum cmake version
cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/main.cpp
    ${PROJECT_SOURCE_DIR}/src/write_tables.cpp
    ${PROJECT_SOURCE_DIR}/src/TwoToTwo_differential_Xsections.cpp
    ${PROJECT_SOURCE_DIR}/src/TwoToTwo_total_Xsections.cpp
  )

add_executable(phoxtrot ${SOURCES})

# execute phoxtrot and generate tables
function(generate_cross_section_tables)
  set(cross_sections "${CMAKE_CURRENT_BINARY_DIR}/tables")
  file(MAKE_DIRECTORY ${cross_sections})

  add_custom_command(
              OUTPUT "${cross_sections}/Total_Xsections_PiRho_PiGamma"
                     "${cross_sections}/Diff_Xsections_PiRho_PiGamma"
                     "${cross_sections}/Total_Xsections_PiPi_RhoGamma"
                     "${cross_sections}/Diff_Xsections_PiPi_RhoGamma"
              COMMAND "${CMAKE_BINARY_DIR}/phoxtrot"
                      "${cross_sections}"
              DEPENDS "${CMAKE_BINARY_DIR}/phoxtrot"
              COMMENT "Calculating cross sections for photon processes."
              VERBATIM
          )
  set(cross_section_files "${cross_sections}/Total_Xsections_PiRho_PiGamma"
                          "${cross_sections}/Diff_Xsections_PiRho_PiGamma"
                          "${cross_sections}/Total_Xsections_PiPi_RhoGamma"
                          "${cross_sections}/Diff_Xsections_PiPi_RhoGamma"
                          PARENT_SCOPE)
  set(Tables_dir "${cross_sections}" PARENT_SCOPE)
endfunction()


# Plot cross sections from phoxtrot output
function(plot_cross_sections)

  set(plots_dir "${CMAKE_CURRENT_BINARY_DIR}/plots")
  file(MAKE_DIRECTORY ${plots_dir})

  add_custom_command(
              OUTPUT "${plots_dir}/Total_Xsec_PiPi_RhoGamma.pdf"
                     "${plots_dir}/Total_Xsec_PiRho_PiGamma.pdf"
                     "${plots_dir}/Diff_Xsec_PiPi_RhoGamma.pdf"
                     "${plots_dir}/Diff_Xsec_PiRho_PiGamma.pdf"
              COMMAND "python2" "${PROJECT_SOURCE_DIR}/src/plot_cross_sections.py"
                      "-d" "${Tables_dir}"
                      "-o" ${plots_dir}
              DEPENDS "${PROJECT_SOURCE_DIR}/src/plot_cross_sections.py"
                      "${cross_section_files}"
              COMMENT "Plotting cross sections for photon processes."
              VERBATIM
          )

    set(cross_section_plots "${plots_dir}/Total_Xsec_PiPi_RhoGamma.pdf"
                            "${plots_dir}/Total_Xsec_PiRho_PiGamma.pdf"
                            "${plots_dir}/Diff_Xsec_PiPi_RhoGamma.pdf"
                            "${plots_dir}/Diff_Xsec_PiRho_PiGamma.pdf"
                            PARENT_SCOPE)
endfunction()


generate_cross_section_tables()
plot_cross_sections()

add_custom_target(plots
    DEPENDS ${cross_section_plots}
)
