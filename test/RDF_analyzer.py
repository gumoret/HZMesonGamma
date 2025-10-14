import ROOT
import argparse

#Following bools are given as input
verbose       = False
isPhiAnalysis = False # for H -> Phi Gamma
isRhoAnalysis = False # for H -> Rho Gamma

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description="RDataFrame analyzer for H→meson+γ")
p.add_argument("meson_option", help="Type <<rho>> for rho, <<phi>> for phi")
p.add_argument("runningOnData_option", help="Type <<signal>> for signal, <<data>> for data")
p.add_argument("rootfile_name", help="Input nanoAOD ROOT file")
p.add_argument("outputfile_option", help="Output ROOT file")
args = p.parse_args()

if args.meson_option == "phi": isPhiAnalysis = True
elif args.meson_option == "rho": isRhoAnalysis = True
else: print("meson_option must be <<phi>> or <<rho>> or <<K*>> or <<DO*>>")


if args.runningOnData_option == "signal": runningOnData = False
elif args.runningOnData_option == "data": runningOnData = True
else: print("runninOnData must be <<signal>> or <<data>>")


input_file = args.rootfile_name
output_file = args.outputfile_option
input_tree_name = "Events"

# Create the dataframe
df = ROOT.RDataFrame(input_tree_name, input_file)

# ---------------------------------------------------------------------
# Trigger
# ---------------------------------------------------------------------
df_trigger = df.Filter("HLT_Photon35_TwoProngs35", "Two-prong photon trigger")

# ---------------------------------------------------------------------
# Meson mass 
# ---------------------------------------------------------------------
if isPhiAnalysis:
    df_meson = df_trigger.Define("mesonTreeMass", "phi_kin_mass[0]")  
elif isRhoAnalysis:
    df_meson = df_trigger.Define("mesonTreeMass", "rho_kin_mass[0]")


# ---------------------------------------------------------------------
# output saving
# ---------------------------------------------------------------------
columns_to_save = ["HLT_Photon35_TwoProngs35", "mesonTreeMass"]

df_meson.Snapshot("tree_output", output_file, columns_to_save)

print(f"Output saved in: {output_file}")
print(f"Tree name: tree_output")