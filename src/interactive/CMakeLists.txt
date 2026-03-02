# TODO: Update to sfml v3.0
# add_executable(SFMLTest main_sfml.cc Renderer.cc)
# target_link_libraries(SFMLTest PRIVATE project_options
#                                PUBLIC ImGui-SFML::ImGui-SFML geom tbb)

add_executable(BeamViz main_beamviz.cc)
target_link_libraries(BeamViz PRIVATE project_options
                              PUBLIC ImGui-SFML::ImGui-SFML geom OpenGL::GL)

add_executable(GeoDisplay main_geodisplay.cc)
target_link_libraries(GeoDisplay PRIVATE project_options
                                 PUBLIC ImGui-SFML::ImGui-SFML geom OpenGL::GL)
