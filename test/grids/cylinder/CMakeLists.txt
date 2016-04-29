# Copy block grid files used for testing
file(GLOB GRIDS "*.h5")
file(COPY ${GRIDS} DESTINATION ${CMAKE_BINARY_DIR}/bin)


