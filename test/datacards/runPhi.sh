combineCards.py GFcat_bdt0=datacard_STAT_Phi_GFcat_bdt0_2024.txt GFcat_bdt1=datacard_STAT_Phi_GFcat_bdt1_2024.txt > datacard_STAT_Phi_GFcat_Comb_2024.txt
echo "### Limits" > ./UL_Phi_GFcat.log
combine datacard_STAT_Phi_GFcat_Comb_2024.txt -M AsymptoticLimits -m 125 -t -1 --run expected --name PhiGFcat >> ./UL_Phi_GFcat.log
echo "### Background-only fit checks - Asinov" >> ./UL_Phi_GFcat.log
combine datacard_STAT_Phi_GFcat_Comb_2024.txt -M FitDiagnostics -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 0 --name PhiGFcat_bkg_only >> ./UL_Phi_GFcat.log
echo "### Total fit checks - Asinov" >> ./UL_Phi_GFcat.log
combine datacard_STAT_Phi_GFcat_Comb_2024.txt -M FitDiagnostics -m 125 --rMin -10 --rMax 10 -t -1 --expectSignal 1 --name PhiGFcat_total >> ./UL_Phi_GFcat.log
echo "### Total fit checks - Toys" >> ./UL_Phi_GFcat.log
combine datacard_STAT_Phi_GFcat_Comb_2024.txt -M FitDiagnostics -m 125 --rMin -10 --rMax 10 -t 500 --expectSignal 1 --name PhiGFcat_toys >> ./UL_Phi_GFcat.log
echo "### Plot impacts" >> ./UL_Phi_GFcat.log
text2workspace.py -m 125 datacard_STAT_Phi_GFcat_Comb_2024.txt
combineTool.py -M Impacts -d datacard_STAT_Phi_GFcat_Comb_2024.root -m 125 --rMin -10 --rMax 10 --doInitialFit --robustFit 1 >> ./UL_Phi_GFcat.log
combineTool.py -M Impacts -d datacard_STAT_Phi_GFcat_Comb_2024.root -m 125 --rMin -10 --rMax 10 --doFits --robustFit 1 >> ./UL_Phi_GFcat.log
combineTool.py -M Impacts -d datacard_STAT_Phi_GFcat_Comb_2024.root -m 125 -o impacts.json
plotImpacts.py -i impacts.json -o impacts_Phi
cp impacts_Phi.pdf /eos/user/c/covarell/www/AN-26-106/
echo "### Plot toys" >> ./UL_Phi_GFcat.log
root -l runPhiToys.C
cp rpull.pdf /eos/user/c/covarell/www/AN-26-106/rpull_Phi.pdf
cp llrdist.pdf /eos/user/c/covarell/www/AN-26-106/llrdist_Phi.pdf
