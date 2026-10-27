combineCards.py GFcat_bdt0=datacard_STAT_Dstar_GFcat_bdt0_2024.txt GFcat_bdt1=datacard_STAT_Dstar_GFcat_bdt1_2024.txt > datacard_STAT_Dstar_GFcat_Comb_2024.txt
echo "### Limits" > ./UL_Dstar_GFcat.log
combine datacard_STAT_Dstar_GFcat_Comb_2024.txt -M AsymptoticLimits -m 125 -t -1 --run expected --name DstarGFcat >> ./UL_Dstar_GFcat.log
echo "### Background-only fit checks - Asinov" >> ./UL_Dstar_GFcat.log
combine datacard_STAT_Dstar_GFcat_Comb_2024.txt -M FitDiagnostics -m 125 --rMin -10 -t -1 --expectSignal 0 --name DstarGFcat_bkg_only >> ./UL_Dstar_GFcat.log
echo "### Total fit checks - Asinov" >> ./UL_Dstar_GFcat.log
combine datacard_STAT_Dstar_GFcat_Comb_2024.txt -M FitDiagnostics -m 125 --rMin -10 -t -1 --expectSignal 1 --name DstarGFcat_total >> ./UL_Dstar_GFcat.log
echo "### Total fit checks - Toys" >> ./UL_Dstar_GFcat.log
combine datacard_STAT_Dstar_GFcat_Comb_2024.txt -M FitDiagnostics -m 125 --rMin -3 -t 500 --expectSignal 1 --name DstarGFcat_toys >> ./UL_Dstar_GFcat.log
echo "### Plot impacts" >> ./UL_Dstar_GFcat.log
text2workspace.py -m 125 datacard_STAT_Dstar_GFcat_Comb_2024.txt
combineTool.py -M Impacts -d datacard_STAT_Dstar_GFcat_Comb_2024.root -m 125 --rMin -10 --doInitialFit --robustFit 1 >> ./UL_Dstar_GFcat.log
combineTool.py -M Impacts -d datacard_STAT_Dstar_GFcat_Comb_2024.root -m 125 --rMin -10 --doFits --robustFit 1 >> ./UL_Dstar_GFcat.log
combineTool.py -M Impacts -d datacard_STAT_Dstar_GFcat_Comb_2024.root -m 125 -o impacts.json
plotImpacts.py -i impacts.json -o impacts
cp impacts.pdf /eos/user/c/covarell/www/AN-26-106/
echo "### Plot toys" >> ./UL_Dstar_GFcat.log
root -l runDstarToys.C
cp rpull.pdf /eos/user/c/covarell/www/AN-26-106/
cp llrdist.pdf /eos/user/c/covarell/www/AN-26-106/
