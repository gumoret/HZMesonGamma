#!/bin/bash
python3 HZMesonGamma.py phi signal rootfiles/input_files/ZPhiGamma.root rootfiles/latest_productions/Signal/ZPhiGamma_output.root
python3 HZMesonGamma.py rho signal rootfiles/input_files/ZRhoGamma.root rootfiles/latest_productions/Signal/ZRhoGamma_output.root
python3 HZMesonGamma.py K signal rootfiles/input_files/ZKstGamma.root rootfiles/latest_productions/Signal/ZKstGamma_output.root
######################################################
python3 HZMesonGamma.py phi signal rootfiles/input_files/HPhiGamma.root rootfiles/latest_productions/Signal/HPhiGamma_output.root
python3 HZMesonGamma.py rho signal rootfiles/input_files/HRhoGamma.root rootfiles/latest_productions/Signal/HRhoGamma_output.root
python3 HZMesonGamma.py K signal rootfiles/input_files/HKstGamma.root rootfiles/latest_productions/Signal/HKstGamma_output.root
#######################################################
######################################################
python3 generate_histos.py Z phi rootfiles/latest_productions/Signal/ZPhiGamma_output.root histos/Signal/ZPhiGamma_signal.root
python3 generate_histos.py Z rho rootfiles/latest_productions/Signal/ZRhoGamma_output.root histos/Signal/ZRhoGamma_signal.root
python3 generate_histos.py Z K rootfiles/latest_productions/Signal/ZKstGamma_output.root histos/Signal/ZKstGamma_signal.root
#####################################################
python3 generate_histos.py H phi rootfiles/latest_productions/Signal/HPhiGamma_output.root histos/Signal/HPhiGamma_signal.root
python3 generate_histos.py H rho rootfiles/latest_productions/Signal/HRhoGamma_output.root histos/Signal/HRhoGamma_signal.root
python3 generate_histos.py H K rootfiles/latest_productions/Signal/HKstGamma_output.root histos/Signal/HKstGamma_signal.root
#######################################################
#######################################################
python3 plot.py Z phi signal SR histos/Signal/ZPhiGamma_signal.root
python3 plot.py Z rho signal SR histos/Signal/ZRhoGamma_signal.root
python3 plot.py Z K signal SR histos/Signal/ZKstGamma_signal.root
######################################################
python3 plot.py H phi signal SR histos/Signal/HPhiGamma_signal.root
python3 plot.py H rho signal SR histos/Signal/HRhoGamma_signal.root
python3 plot.py H K signal SR histos/Signal/HKstGamma_signal.root