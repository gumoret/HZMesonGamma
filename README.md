# HZMesonGamma

1) HZMesonGamma.py
python3 HZMesonGamma.py <<meson>> <<signal/data>> <<nano input file>> <<output rootfile>>
NB: take the nano inputfiles from /eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/Run3-D05/, and do the hadd command to save the rootfiles in a single file
2) generate_histos.py
python3 generate_histos.py <<boson>> <<meson>> <<input rootfile>> <<output histogram rootfile>>
NB: the <<input rootfile>> is the <<output rootfile>> from HZMesonGamma.py
3) plot.py
python3 plot.py <<boson>> <<meson>> <<signal/data>> <<SR/CR>> <<input histogram rootfile>>
NB: the <<input histogram rootfile>> is the <<output histogram rootfile>> from generate_histos.py




