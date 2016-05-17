#!/bin/bash

# Build documentation
doxygen ChiDG.cfg


# Copy supporting files to HTML directory
cp -r extra_files/* html/


