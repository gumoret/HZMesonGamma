#!/bin/bash
python3 TriggerAnalysis.py phi signal rootfiles/input_files/HPhiGamma.root histos/Signal/HPhiGamma_triggers.root TwoProngs
python3 TriggerAnalysis.py phi signal rootfiles/input_files/HPhiGamma.root histos/Signal/HPhiGamma_triggers.root Photon50
python3 TriggerAnalysis.py phi signal rootfiles/input_files/HPhiGamma.root histos/Signal/HPhiGamma_triggers.root OR 
python3 plot_trigger.py H phi histos/Signal/HPhiGamma_triggers.root
#####################################################
python3 TriggerAnalysis.py rho signal rootfiles/input_files/HRhoGamma.root histos/Signal/HRhoGamma_triggers.root TwoProngs
python3 TriggerAnalysis.py rho signal rootfiles/input_files/HRhoGamma.root histos/Signal/HRhoGamma_triggers.root Photon50
python3 TriggerAnalysis.py rho signal rootfiles/input_files/HRhoGamma.root histos/Signal/HRhoGamma_triggers.root OR 
python3 plot_trigger.py H rho histos/Signal/HRhoGamma_triggers.root
#####################################################
python3 TriggerAnalysis.py K signal rootfiles/input_files/HKstGamma.root histos/Signal/HKstGamma_triggers.root TwoProngs
python3 TriggerAnalysis.py K signal rootfiles/input_files/HKstGamma.root histos/Signal/HKstGamma_triggers.root Photon50
python3 TriggerAnalysis.py K signal rootfiles/input_files/HKstGamma.root histos/Signal/HKstGamma_triggers.root OR 
python3 plot_trigger.py H K histos/Signal/HKstGamma_triggers.root
#######################################################
#######################################################
python3 TriggerAnalysis.py phi signal rootfiles/input_files/ZPhiGamma.root histos/Signal/ZPhiGamma_triggers.root TwoProngs
python3 TriggerAnalysis.py phi signal rootfiles/input_files/ZPhiGamma.root histos/Signal/ZPhiGamma_triggers.root Photon50
python3 TriggerAnalysis.py phi signal rootfiles/input_files/ZPhiGamma.root histos/Signal/ZPhiGamma_triggers.root OR 
python3 plot_trigger.py Z phi histos/Signal/ZPhiGamma_triggers.root
#####################################################
python3 TriggerAnalysis.py rho signal rootfiles/input_files/ZRhoGamma.root histos/Signal/ZRhoGamma_triggers.root TwoProngs
python3 TriggerAnalysis.py rho signal rootfiles/input_files/ZRhoGamma.root histos/Signal/ZRhoGamma_triggers.root Photon50
python3 TriggerAnalysis.py rho signal rootfiles/input_files/ZRhoGamma.root histos/Signal/ZRhoGamma_triggers.root OR 
python3 plot_trigger.py Z rho histos/Signal/ZRhoGamma_triggers.root
#####################################################
python3 TriggerAnalysis.py K signal rootfiles/input_files/ZKstGamma.root histos/Signal/ZKstGamma_triggers.root TwoProngs
python3 TriggerAnalysis.py K signal rootfiles/input_files/ZKstGamma.root histos/Signal/ZKstGamma_triggers.root Photon50
python3 TriggerAnalysis.py K signal rootfiles/input_files/ZKstGamma.root histos/Signal/ZKstGamma_triggers.root OR 
python3 plot_trigger.py Z K histos/Signal/ZKstGamma_triggers.root
