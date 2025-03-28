# HZMesonGamma Analysis
  
Below are the instructions for running the main scripts.

---

## 1. HZMesonGamma.py
```bash
python3 HZMesonGamma.py <meson> <signal/data> <nano input file> <output rootfile>
```
**Note:**  
Retrieve the `<nano input file>` from the directory:  
`/eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/Run3-D05/`  
Use the `hadd` command to merge multiple NanoAOD files into a single ROOT file.

---

## 2. generate_histos.py
```bash
python3 generate_histos.py <boson> <meson> <input rootfile> <output histogram rootfile>
```
**Note:**  
The `<input rootfile>` is the `<output rootfile>` generated by `HZMesonGamma.py`.

---

## 3. plot.py
```bash
python3 plot.py <boson> <meson> <signal/data> <SR/CR> <input histogram rootfile>
```
**Note:**  
The `<input histogram rootfile>` is the `<output histogram rootfile>` generated by `generate_histos.py`.

---








