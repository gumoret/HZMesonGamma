#!/bin/bash
python3 generate_histos.py data Z phi CR unblind preselection rootfiles/latest_productions/HPhiGamma2022Tau.root histos/Z/CR_ZPhi2022Tau_preselection_Sidebands.root
python3 generate_histos.py data Z phi SR blind preselection rootfiles/latest_productions/HPhiGamma2022Tau.root histos/Z/SR_ZPhi2022Tau_preselection_Data.root
python3 generate_histos.py signal Z phi SR unblind preselection rootfiles/latest_productions/ZPhiGamma_Signal.root histos/Z/SR_ZPhi_preselection_Signal.root

python3 normalize_for_BkgEstimation.py Z histos/Z/SR_ZPhi2022Tau_preselection_Data.root histos/Z/CR_ZPhi2022Tau_preselection_Sidebands.root 
python3 plot_histos_bkgEstimation.py 0.05 loose Z phi histos/Z/SR_ZPhi2022Tau_preselection_Data.root histos/Z/SR_ZPhi_preselection_Signal.root histos/Z/CR_ZPhi2022Tau_preselection_SidebandsNorm.root 
###
python3 generate_histos.py data Z rho CR unblind preselection rootfiles/latest_productions/HRhoGamma2022Tau.root histos/Z/CR_ZRho2022Tau_preselection_Sidebands.root
python3 generate_histos.py data Z rho SR blind preselection rootfiles/latest_productions/HRhoGamma2022Tau.root histos/Z/SR_ZRho2022Tau_preselection_Data.root
python3 generate_histos.py signal Z rho SR unblind preselection rootfiles/latest_productions/ZRhoGamma_Signal.root histos/Z/SR_ZRho_preselection_Signal.root

python3 normalize_for_BkgEstimation.py Z histos/Z/SR_ZRho2022Tau_preselection_Data.root histos/Z/CR_ZRho2022Tau_preselection_Sidebands.root
python3 plot_histos_bkgEstimation.py 0.05 loose Z rho histos/Z/SR_ZRho2022Tau_preselection_Data.root histos/Z/SR_ZRho_preselection_Signal.root histos/Z/CR_ZRho2022Tau_preselection_SidebandsNorm.root 