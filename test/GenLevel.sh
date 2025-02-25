#!/bin/bash
python3 HZMesonGammaGenLevel.py rootfiles/input_files/HPhiGamma.root rootfiles/latest_productions/HPhiGammaGenLevel_output.root
python3 HZMesonGammaGenLevel.py rootfiles/input_files/HRhoGamma.root rootfiles/latest_productions/HRhoGammaGenLevel_output.root
python3 HZMesonGammaGenLevel.py rootfiles/input_files/HKstGamma.root rootfiles/latest_productions/HKstGammaGenLevel_output.root
###############################################################################################################################
python3 HZMesonGammaGenLevel.py rootfiles/input_files/ZPhiGamma.root rootfiles/latest_productions/ZPhiGammaGenLevel_output.root
python3 HZMesonGammaGenLevel.py rootfiles/input_files/ZRhoGamma.root rootfiles/latest_productions/ZRhoGammaGenLevel_output.root
python3 HZMesonGammaGenLevel.py rootfiles/input_files/ZKstGamma.root rootfiles/latest_productions/ZKstGammaGenLevel_output.root
###############################################################################################################################
###############################################################################################################################
python3 generate_histos_GenLevel.py H phi rootfiles/latest_productions/HPhiGammaGenLevel_output.root histos/GenLevel/HPhiGamma_GenLevel.root
python3 generate_histos_GenLevel.py H rho rootfiles/latest_productions/HRhoGammaGenLevel_output.root histos/GenLevel/HRhoGamma_GenLevel.root
python3 generate_histos_GenLevel.py H K rootfiles/latest_productions/HKstGammaGenLevel_output.root histos/GenLevel/HKstGamma_GenLevel.root
###############################################################################################################################
python3 generate_histos_GenLevel.py Z phi rootfiles/latest_productions/ZPhiGammaGenLevel_output.root histos/GenLevel/ZPhiGamma_GenLevel.root
python3 generate_histos_GenLevel.py Z rho rootfiles/latest_productions/ZRhoGammaGenLevel_output.root histos/GenLevel/ZRhoGamma_GenLevel.root
python3 generate_histos_GenLevel.py Z K rootfiles/latest_productions/ZKstGammaGenLevel_output.root histos/GenLevel/ZKstGamma_GenLevel.root
###############################################################################################################################
###############################################################################################################################
python3 plot_GenLevel.py H phi histos/GenLevel/HPhiGamma_GenLevel.root
python3 plot_GenLevel.py H rho histos/GenLevel/HRhoGamma_GenLevel.root
python3 plot_GenLevel.py H K histos/GenLevel/HKstGamma_GenLevel.root
###############################################################################################################################
python3 plot_GenLevel.py Z phi histos/GenLevel/ZPhiGamma_GenLevel.root
python3 plot_GenLevel.py Z rho histos/GenLevel/ZRhoGamma_GenLevel.root
python3 plot_GenLevel.py Z K histos/GenLevel/ZKstGamma_GenLevel.root