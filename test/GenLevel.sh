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
python3 generate_histos_GenLevel.py H phi rootfiles/latest_productions/HPhiGammaGenLevel_output.root histos/HPhiGamma_GenLevel.root
python3 generate_histos_GenLevel.py H rho rootfiles/latest_productions/HRhoGammaGenLevel_output.root histos/HRhoGamma_GenLevel.root
python3 generate_histos_GenLevel.py H K rootfiles/latest_productions/HKstGammaGenLevel_output.root histos/HKstGamma_GenLevel.root
###############################################################################################################################
python3 generate_histos_GenLevel.py Z phi rootfiles/latest_productions/ZPhiGammaGenLevel_output.root histos/ZPhiGamma_GenLevel.root
python3 generate_histos_GenLevel.py Z rho rootfiles/latest_productions/ZRhoGammaGenLevel_output.root histos/ZRhoGamma_GenLevel.root
python3 generate_histos_GenLevel.py Z K rootfiles/latest_productions/ZKstGammaGenLevel_output.root histos/ZKstGamma_GenLevel.root
###############################################################################################################################
###############################################################################################################################
python3 plot_GenLevel.py H phi histos/HPhiGamma_GenLevel.root
python3 plot_GenLevel.py H rho histos/HRhoGamma_GenLevel.root
python3 plot_GenLevel.py H K histos/HKstGamma_GenLevel.root
###############################################################################################################################
python3 plot_GenLevel.py Z phi histos/ZPhiGamma_GenLevel.root
python3 plot_GenLevel.py Z rho histos/ZRhoGamma_GenLevel.root
python3 plot_GenLevel.py Z K histos/ZKstGamma_GenLevel.root