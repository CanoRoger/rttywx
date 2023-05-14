add_library(librttywx SHARED fftfilt.cxx rtty_rx.cpp)

add_executable(rtty_rx_from_file rtty_rx_from_file.cpp)
target_link_libraries(rtty_rx_from_file librttywx)

include(GNUInstallDirs)
install(TARGETS librttywx rtty_rx_from_file)
install(FILES rtty_rx.h fftfilt.h gfft.h complex.h TYPE INCLUDE)
