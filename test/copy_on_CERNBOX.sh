#!/bin/bash

# Controllo argomenti
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <INPUT_PATH> <OUTPUT_PATH>"
    exit 1
fi

INPUT_PATH="$1"
OUTPUT_PATH="$2"

LXPLUS_HOST="eferrand@lxplus.cern.ch:/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/rootfiles/latest_productions"

echo "Copying:"
echo "  From: $INPUT_PATH"
echo "  To:   $LXPLUS_HOST/$OUTPUT_PATH"
echo ""

rsync -avz "$INPUT_PATH" "$LXPLUS_HOST/$OUTPUT_PATH"