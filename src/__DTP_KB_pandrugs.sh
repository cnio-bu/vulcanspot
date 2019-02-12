#!/bin/bash

### TITLE : Knokledge-based Drug prescription using PanDrugs
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es

# Prepare input for PanDrugs prescription
python ./src/PanDrugs/input_file.py
# Drug Prescription
python ./src/PanDrugs/process_drugs.py
