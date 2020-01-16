#!/bin/bash

INPUT_DIRECTORY=unfolding/input/

FILE="${INPUT_DIRECTORY}unfold_spectrum.cfg"
if [ ! -f "$FILE" ]; then
    cp "${INPUT_DIRECTORY}template_unfold_spectrum.cfg" "$FILE"
fi

FILE="${INPUT_DIRECTORY}unfold_trend.cfg"
if [ ! -f "$FILE" ]; then
    cp "${INPUT_DIRECTORY}template_unfold_trend.cfg" "$FILE"
fi

FILE="${INPUT_DIRECTORY}plot_spectra.cfg"
if [ ! -f "$FILE" ]; then
    cp "${INPUT_DIRECTORY}template_plot_spectra.cfg" "$FILE"
fi

FILE="${INPUT_DIRECTORY}plot_lines.cfg"
if [ ! -f "$FILE" ]; then
    cp "${INPUT_DIRECTORY}template_plot_lines.cfg" "$FILE"
fi

echo "Initialization complete"