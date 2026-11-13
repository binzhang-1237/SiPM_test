cmake_minimum_required(VERSION 3.16)
project(SiPMTest LANGUAGES CXX)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad Spectrum Minuit)

add_executable(analyze_sipm analyze_sipm.cpp)
target_compile_features(analyze_sipm PRIVATE cxx_std_17)
target_compile_options(analyze_sipm PRIVATE -Wall -Wextra)
target_link_libraries(analyze_sipm PRIVATE ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Spectrum ROOT::Minuit)
