#!/bin/bash
python3 generate_histos.py data H phi CR unblind preselection rootfiles/latest_productions/HPhiGamma2022Tau.root histos/H/CR_HPhi2022Tau_preselection_Sidebands.root
python3 generate_histos.py data H phi SR blind preselection rootfiles/latest_productions/HPhiGamma2022Tau.root histos/H/SR_HPhi2022Tau_preselection_Data.root
python3 generate_histos.py signal H phi SR unblind preselection rootfiles/latest_productions/HPhiGamma_Signal.root histos/H/SR_HPhi_preselection_Signal.root

python3 normalize_for_BkgEstimation.py H histos/H/SR_HPhi2022Tau_preselection_Data.root histos/H/CR_HPhi2022Tau_preselection_Sidebands.root 
python3 plot_histos_bkgEstimation.py 0.05 loose H phi histos/H/SR_HPhi2022Tau_preselection_Data.root histos/H/SR_HPhi_preselection_Signal.root histos/H/CR_HPhi2022Tau_preselection_SidebandsNorm.root 
###
python3 generate_histos.py data H rho CR unblind preselection rootfiles/latest_productions/HRhoGamma2022Tau.root histos/H/CR_HRho2022Tau_preselection_Sidebands.root
python3 generate_histos.py data H rho SR blind preselection rootfiles/latest_productions/HRhoGamma2022Tau.root histos/H/SR_HRho2022Tau_preselection_Data.root
python3 generate_histos.py signal H rho SR unblind preselection rootfiles/latest_productions/HRhoGamma_Signal.root histos/H/SR_HRho_preselection_Signal.root

python3 normalize_for_BkgEstimation.py H histos/H/SR_HRho2022Tau_preselection_Data.root histos/H/CR_HRho2022Tau_preselection_Sidebands.root
python3 plot_histos_bkgEstimation.py 0.05 loose H rho histos/H/SR_HRho2022Tau_preselection_Data.root histos/H/SR_HRho_preselection_Signal.root histos/H/CR_HRho2022Tau_preselection_SidebandsNorm.root 
###
python3 generate_histos.py data H K CR unblind preselection rootfiles/latest_productions/HKstGamma2022Tau.root histos/H/CR_HKst2022Tau_preselection_Sidebands.root
python3 generate_histos.py data H K SR blind preselection rootfiles/latest_productions/HKstGamma2022Tau.root histos/H/SR_HKst2022Tau_preselection_Data.root
python3 generate_histos.py signal H K SR unblind preselection rootfiles/latest_productions/HKstGamma_Signal.root histos/H/SR_HKst_preselection_Signal.root

python3 normalize_for_BkgEstimation.py H histos/H/SR_HKst2022Tau_preselection_Data.root histos/H/CR_HKst2022Tau_preselection_Sidebands.root
python3 plot_histos_bkgEstimation.py 0.05 loose H K histos/H/SR_HKst2022Tau_preselection_Data.root histos/H/SR_HKst_preselection_Signal.root histos/H/CR_HKst2022Tau_preselection_SidebandsNorm.root 