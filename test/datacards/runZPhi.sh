echo "### Limits" > ./UL_Phi_Z.log
combine datacard_STAT_Phi_Z_bdt0_2024.txt -M AsymptoticLimits -m 125 -t -1 --run expected --name PhiZ >> ./UL_Phi_Z.log
echo "### Background-only fit checks - Asinov" >> ./UL_Phi_Z.log
combine datacard_STAT_Phi_Z_bdt0_2024.txt -M FitDiagnostics -m 125 --rMin -10 -t -1 --expectSignal 0 --name PhiZ_bkg_only >> ./UL_Phi_Z.log
echo "### Total fit checks - Asinov" >> ./UL_Phi_Z.log
combine datacard_STAT_Phi_Z_bdt0_2024.txt -M FitDiagnostics -m 125 --rMin -10 -t -1 --expectSignal 1 --name PhiZ_total >> ./UL_Phi_Z.log
echo "### Total fit checks - Toys" >> ./UL_Phi_Z.log
combine datacard_STAT_Phi_Z_bdt0_2024.txt -M FitDiagnostics -m 125 --rMin -3 -t 500 --expectSignal 1 --name PhiZ_toys >> ./UL_Phi_Z.log
echo "### Plot impacts" >> ./UL_Phi_Z.log
text2workspace.py -m 125 datacard_STAT_Phi_Z_bdt0_2024.txt
combineTool.py -M Impacts -d datacard_STAT_Phi_Z_bdt0_2024.root -m 125 --rMin -10 --doInitialFit --robustFit 1 >> ./UL_Phi_Z.log
combineTool.py -M Impacts -d datacard_STAT_Phi_Z_bdt0_2024.root -m 125 --rMin -10 --doFits --robustFit 1 >> ./UL_Phi_Z.log
combineTool.py -M Impacts -d datacard_STAT_Phi_Z_bdt0_2024.root -m 125 -o impacts.json
plotImpacts.py -i impacts.json -o impacts_ZPhi
cp impacts_ZPhi.pdf /eos/user/c/covarell/www/AN-26-106/
echo "### Plot toys" >> ./UL_Phi_Z.log
root -l runZPhiToys.C
cp rpull.pdf /eos/user/c/covarell/www/AN-26-106/rpull_ZPhi.pdf
cp llrdist.pdf /eos/user/c/covarell/www/AN-26-106/llrdist_ZPhi.pdf
