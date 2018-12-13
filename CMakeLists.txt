# project name
project(phoxtrot)

# minimum cmake version
cmake_minimum_required(VERSION 2.8.11)

include_directories("${PROJECT_SOURCE_DIR}/src")

# compile code and create executable 'phoxtrot'
add_executable(phoxtrot src/main.cpp)

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
endfunction()


generate_cross_section_tables()
add_custom_target(phoxtrot_tables
    DEPENDS ${cross_section_files}
)
