"""Standard InChIs shared by the tests, named after the structure they describe.

Where the InChI Technical Manual prints a standard InChI that fits a case, it is used as printed and the figure is
cited: S. E. Stein, S. R. Heller, D. V. Tchekhovskoi and I. V. Pletnev, "IUPAC International Chemical Identifier
(InChI), InChI version 1, software version 1.04 (2011), Technical Manual". Page numbers are PDF pages. The other
InChIs are the standard InChIs RDKit generates from the SMILES given next to them.
"""

# --- From the InChI Technical Manual ---------------------------------------------------------------------------------

# Fig. 25: L-(+)-tartaric acid (2R,3R), and its enantiomer, which InChI writes with the inversion flag /m0 instead of /m1
TARTARIC_L = "InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/t1-,2-/m1/s1"
TARTARIC_D = "InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/t1-,2-/m0/s1"
# PDF p. 50: relative stereo (/s2) only exists in non-standard InChI (prefix "InChI=1/"); standard InChI allows "only
# absolute stereo or no stereo at all" (PDF p. 15)
TARTARIC_RELATIVE = "InChI=1/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/t1-,2-/s2"
# Fig. 20: butane-2,3-diol with both centres defined (a), and with centre 4 undefined, written "4?" (b)
BUTANEDIOL_TWO_CENTRES = "InChI=1S/C4H10O2/c1-3(5)4(2)6/h3-6H,1-2H3/t3-,4-/m0/s1"
BUTANEDIOL_ONE_CENTRE = "InChI=1S/C4H10O2/c1-3(5)4(2)6/h3-6H,1-2H3/t3-,4?/m0/s1"
# Fig. 28: a triene whose structures (a) and (b) differ in the geometry of one of three double bonds
TRIENE_A = "InChI=1S/C9H12/c1-4-7-8(5-2)9(7)6-3/h4-6H,1-3H3/b7-4-,8-5-,9-6-"
TRIENE_B = "InChI=1S/C9H12/c1-4-7-8(5-2)9(7)6-3/h4-6H,1-3H3/b7-4-,8-5+,9-6-"
# PDF p. 65: trans-but-2-ene
TRANS_BUT_2_ENE = "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+"
# Fig. 32: glycine; its neutral and zwitterionic forms share one standard InChI, the anion and cation add /p
GLYCINE = "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)"
GLYCINATE = "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)/p-1"
GLYCINIUM = "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)/p+1"
# Fig. 16 (PDF p. 39): a cation drawn with the charge in two places (16c, 16f), which the manual says "should have same
# identifiers. Unfortunately, they don't" (their mobile-H groups differ), and the neutral 16g
CATION_16C = "InChI=1S/C9H9N3OS/c1-12(2)9-10-5-3-7(13)8(14)4-6(5)11-9/h3-4H,1-2H3,(H,10,11,13,14)/p+1"
CATION_16F = "InChI=1S/C9H9N3OS/c1-12(2)9-10-5-3-7(13)8(14)4-6(5)11-9/h3-4H,1-2H3,(H,10,11)/p+1"
NEUTRAL_16G = "InChI=1S/C9H9N3OS/c1-12(2)9-10-5-3-7(13)8(14)4-6(5)11-9/h3-4H,1-2H3,(H,10,11)"
# Fig. 17: benzene with two 13C atoms and one D; unlabelled benzene is the same InChI without the /i layer
BENZENE_LABELLED = "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H/i1+1,4+1D"
BENZENE = "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"
# Fig. 19: urea with two exchangeable D, which are written only in the isotopic /h sublayer ("/i/hD2")
UREA_D2 = "InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)/i/hD2"
UREA = "InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)"
# Fig. 31: caffeine
CAFFEINE = "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"

# --- Generated with RDKit ----------------------------------------------------------------------------------------------

MESO_TARTARIC = "InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/t1-,2+"  # O=C(O)[C@H](O)[C@H](O)C(=O)O
# threonine with one of its two centres defined: C2 (bearing the OH) or C3 (bearing the NH2), and C3 protonated
THREONINE_C2_ONLY = "InChI=1S/C4H9NO3/c1-2(6)3(5)4(7)8/h2-3,6H,5H2,1H3,(H,7,8)/t2-,3?/m1/s1"  # C[C@@H](O)C(N)C(=O)O
THREONINE_C3_ONLY = "InChI=1S/C4H9NO3/c1-2(6)3(5)4(7)8/h2-3,6H,5H2,1H3,(H,7,8)/t2?,3-/m0/s1"  # CC(O)[C@H](N)C(=O)O
THREONINIUM_C3_ONLY = (  # CC(O)[C@H]([NH3+])C(=O)O
    "InChI=1S/C4H9NO3/c1-2(6)3(5)4(7)8/h2-3,6H,5H2,1H3,(H,7,8)/p+1/t2?,3-/m0/s1"
)
# C3 only, with the opposite configuration to L-threonine's: D-threonine or L-allothreonine
THREONINE_C3_ONLY_INVERTED = "InChI=1S/C4H9NO3/c1-2(6)3(5)4(7)8/h2-3,6H,5H2,1H3,(H,7,8)/t2?,3-/m1/s1"  # CC(O)[C@@H](N)C(=O)O
# hex-4-ene-2,3-diol with only its double bond defined, only its two centres, or all three
HEXENEDIOL_E_ONLY = "InChI=1S/C6H12O2/c1-3-4-6(8)5(2)7/h3-8H,1-2H3/b4-3+"  # C/C=C/C(O)C(C)O
HEXENEDIOL_CENTRES_ONLY = (  # CC=C[C@@H](O)[C@H](C)O
    "InChI=1S/C6H12O2/c1-3-4-6(8)5(2)7/h3-8H,1-2H3/t5-,6+/m0/s1"
)
HEXENEDIOL_CENTRES_ONLY_INVERTED = (  # CC=C[C@H](O)[C@@H](C)O
    "InChI=1S/C6H12O2/c1-3-4-6(8)5(2)7/h3-8H,1-2H3/t5-,6+/m1/s1"
)
HEXENEDIOL_E_AND_CENTRES = (  # C/C=C/[C@@H](O)[C@H](C)O
    "InChI=1S/C6H12O2/c1-3-4-6(8)5(2)7/h3-8H,1-2H3/b4-3+/t5-,6+/m0/s1"
)
# hexa-2,4-diene with one of its double bonds defined, which lists the other as undefined ("6-4?")
HEXADIENE_ONE_DOUBLE_BOND = "InChI=1S/C6H10/c1-3-5-6-4-2/h3-6H,1-2H3/b5-3+,6-4?"  # C/C=C/C=CC
# 2,3,4-trihydroxyglutaric acid with C3 undefined: C3 is stereogenic only when C2 and C4 differ, (2R,4S), so only
# that form lists it ("1?")
TRIHYDROXYGLUTARIC_2R_4R = (  # OC(=O)[C@H](O)C(O)[C@@H](O)C(=O)O
    "InChI=1S/C5H8O7/c6-1(2(7)4(9)10)3(8)5(11)12/h1-3,6-8H,(H,9,10)(H,11,12)/t2-,3-/m1/s1"
)
TRIHYDROXYGLUTARIC_2R_4S = (  # OC(=O)[C@H](O)C(O)[C@H](O)C(=O)O
    "InChI=1S/C5H8O7/c6-1(2(7)4(9)10)3(8)5(11)12/h1-3,6-8H,(H,9,10)(H,11,12)/t1?,2-,3+"
)
# 1-phenylethylamine with L-tartaric acid; with the amine undefined, InChI lists nothing for that component ("t;")
PHENYLETHYLAMINE_R_TARTRATE = (  # C[C@@H](N)c1ccccc1.OC(=O)[C@H](O)[C@@H](O)C(=O)O
    "InChI=1S/C8H11N.C4H6O6/c1-7(9)8-5-3-2-4-6-8;5-1(3(7)8)2(6)4(9)10/h2-7H,9H2,1H3;1-2,5-6H,(H,7,8)(H,9,10)"
    "/t7-;1-,2-/m11/s1"
)
PHENYLETHYLAMINE_TARTRATE = (  # CC(N)c1ccccc1.OC(=O)[C@H](O)[C@@H](O)C(=O)O
    "InChI=1S/C8H11N.C4H6O6/c1-7(9)8-5-3-2-4-6-8;5-1(3(7)8)2(6)4(9)10/h2-7H,9H2,1H3;1-2,5-6H,(H,7,8)(H,9,10)"
    "/t;1-,2-/m.1/s1"
)
L_ALANINE = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"  # C[C@H](N)C(=O)O
D_ALANINE = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1"  # C[C@@H](N)C(=O)O
ALANINE = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)"  # CC(N)C(=O)O
D_ALANINATE = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/p-1/t2-/m1/s1"  # C[C@@H](N)C(=O)[O-]
L_ALANINE_D3 = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1/i1D3"  # [2H]C([2H])([2H])[C@H](N)C(=O)O
L_ALANINE_1_13C = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1/i3+1"  # C[C@H](N)[13C](=O)O
ALANINE_HCL = "InChI=1S/C3H7NO2.ClH/c1-2(4)3(5)6;/h2H,4H2,1H3,(H,5,6);1H"  # CC(N)C(=O)O.Cl
L_ALANINE_HCL = "InChI=1S/C3H7NO2.ClH/c1-2(4)3(5)6;/h2H,4H2,1H3,(H,5,6);1H/t2-;/m0./s1"  # C[C@H](N)C(=O)O.Cl
D_ALANINE_HCL = "InChI=1S/C3H7NO2.ClH/c1-2(4)3(5)6;/h2H,4H2,1H3,(H,5,6);1H/t2-;/m1./s1"  # C[C@@H](N)C(=O)O.Cl
D_ALANINE_HBR = "InChI=1S/C3H7NO2.BrH/c1-2(4)3(5)6;/h2H,4H2,1H3,(H,5,6);1H/t2-;/m1./s1"  # C[C@@H](N)C(=O)O.Br
# two alanines written once with the n* repeat notation (2C3H7NO2, c2*, t2*2-) and one /m flag each (/m00.)
L_ALANINE2_SULFATE = (  # C[C@H](N)C(=O)O.C[C@H](N)C(=O)O.OS(=O)(=O)O
    "InChI=1S/2C3H7NO2.H2O4S/c2*1-2(4)3(5)6;1-5(2,3)4/h2*2H,4H2,1H3,(H,5,6);(H2,1,2,3,4)/t2*2-;/m00./s1"
)
D_ALANINE2_SULFATE = (  # C[C@@H](N)C(=O)O.C[C@@H](N)C(=O)O.OS(=O)(=O)O
    "InChI=1S/2C3H7NO2.H2O4S/c2*1-2(4)3(5)6;1-5(2,3)4/h2*2H,4H2,1H3,(H,5,6);(H2,1,2,3,4)/t2*2-;/m11./s1"
)
# levalbuterol tartrate, (R)-salbutamol x2 with L-tartaric acid, and with D-tartaric acid: the tartrate's /m flag is the
# third (/m001 vs /m000), so it is read correctly only once the n* repeat is expanded
LEVALBUTEROL_L_TARTRATE = (  # CC(C)(C)NC[C@H](O)c1ccc(O)c(CO)c1 x2 . O=C(O)[C@H](O)[C@@H](O)C(=O)O
    "InChI=1S/2C13H21NO3.C4H6O6/c2*1-13(2,3)14-7-12(17)9-4-5-11(16)10(6-9)8-15;5-1(3(7)8)2(6)4(9)10"
    "/h2*4-6,12,14-17H,7-8H2,1-3H3;1-2,5-6H,(H,7,8)(H,9,10)/t2*12-;1-,2-/m001/s1"
)
LEVALBUTEROL_D_TARTRATE = (  # CC(C)(C)NC[C@H](O)c1ccc(O)c(CO)c1 x2 . O=C(O)[C@@H](O)[C@H](O)C(=O)O
    "InChI=1S/2C13H21NO3.C4H6O6/c2*1-13(2,3)14-7-12(17)9-4-5-11(16)10(6-9)8-15;5-1(3(7)8)2(6)4(9)10"
    "/h2*4-6,12,14-17H,7-8H2,1-3H3;1-2,5-6H,(H,7,8)(H,9,10)/t2*12-;1-,2-/m000/s1"
)
# an achiral first component, so the /m layer starts with "." (no flag)
TEA_D_LACTATE = (  # CCN(CC)CC.C[C@@H](O)C(=O)O
    "InChI=1S/C6H15N.C3H6O3/c1-4-7(5-2)6-3;1-2(4)3(5)6/h4-6H2,1-3H3;2,4H,1H3,(H,5,6)/t;2-/m.1/s1"
)
GLYCINE_HCL = "InChI=1S/C2H5NO2.ClH/c3-1-2(4)5;/h1,3H2,(H,4,5);1H"  # NCC(=O)O.Cl
GLYCINE_HBR = "InChI=1S/C2H5NO2.BrH/c3-1-2(4)5;/h1,3H2,(H,4,5);1H"  # NCC(=O)O.Br
ACETIC_ACID = "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"  # CC(=O)O
SODIUM_ACETATE = "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"  # CC(=O)[O-].[Na+]
ACETONE = "InChI=1S/C3H6O/c1-3(2)4/h1-2H3"  # CC(C)=O
PROPEN_2_OL = "InChI=1S/C3H6O/c1-3(2)4/h4H,1H2,2H3"  # C=C(C)O
S_HYDROXYPENTANONE = "InChI=1S/C5H10O2/c1-4(6)3-5(2)7/h4,6H,3H2,1-2H3/t4-/m0/s1"  # CC(=O)C[C@H](C)O
S_PENTENEDIOL = "InChI=1S/C5H10O2/c1-4(6)3-5(2)7/h3-4,6-7H,1-2H3/t4-/m0/s1"  # CC(O)=C[C@H](C)O, its enol
R_PENTENEDIOL = "InChI=1S/C5H10O2/c1-4(6)3-5(2)7/h3-4,6-7H,1-2H3/t4-/m1/s1"  # CC(O)=C[C@@H](C)O
BUT_1_ENE = "InChI=1S/C4H8/c1-3-4-2/h3H,1,4H2,2H3"  # C=CCC
BUT_2_ENE = "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3"  # CC=CC
BUTANE = "InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3"  # CCCC
ETHANOL = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"  # CCO
DIMETHYL_ETHER = "InChI=1S/C2H6O/c1-3-2/h1-2H3"  # COC
O_CHLOROTOLUENE = "InChI=1S/C7H7Cl/c1-6-4-2-3-5-7(6)8/h2-5H,1H3"  # Cc1ccccc1Cl
P_CHLOROTOLUENE = "InChI=1S/C7H7Cl/c1-6-2-4-7(8)5-3-6/h2-5H,1H3"  # Cc1ccc(Cl)cc1
# a hydroxymethyl or a methoxy group on a benzene ring
BENZYL_ALCOHOL = "InChI=1S/C7H8O/c8-6-7-4-2-1-3-5-7/h1-5,8H,6H2"  # OCc1ccccc1
ANISOLE = "InChI=1S/C7H8O/c1-8-7-5-3-2-4-6-7/h2-6H,1H3"  # COc1ccccc1
ETHANE = "InChI=1S/C2H6/c1-2/h1-2H3"  # CC
ETHYL_RADICAL = "InChI=1S/C2H5/c1-2/h1H2,2H3"  # [CH2]C
IBUPROFEN = (  # CC(C)Cc1ccc(cc1)C(C)C(=O)O
    "InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)"
)
IBUPROFEN_METHYL_ESTER = (  # CC(C)Cc1ccc(cc1)C(C)C(=O)OC
    "InChI=1S/C14H20O2/c1-10(2)9-12-5-7-13(8-6-12)11(3)14(15)16-4/h5-8,10-11H,9H2,1-4H3"
)
THEOPHYLLINE = "InChI=1S/C7H8N4O2/c1-10-5-4(8-3-9-5)6(12)11(2)7(10)13/h3H,1-2H3,(H,8,9)"  # Cn1c(=O)c2[nH]cnc2n(C)c1=O
# the deuterium makes C1 a stereocentre: its stereo is only in the isotopic sublayers (/i2D/t2-/m0/s1)
S_ETHANOL_1_D = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i2D/t2-/m0/s1"  # [2H][C@@H](C)O
R_ETHANOL_1_D = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i2D/t2-/m1/s1"  # [2H][C@H](C)O
ETHANOL_1_D = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i2D"  # [2H]C(C)O
# 13C in one carboxyl makes meso-tartaric acid chiral: the main /t layer keeps the unlabelled meso stereo, and the
# isotopic sublayer adds /m and /s. Labelling the other carboxyl gives the mirror image.
MESO_TARTARIC_13C = (  # O=[13C](O)[C@H](O)[C@H](O)C(=O)O
    "InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/t1-,2+/i3+1/m1/s1"
)
MESO_TARTARIC_13C_OTHER_END = (  # O=C(O)[C@H](O)[C@H](O)[13C](=O)O
    "InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/t1-,2+/i3+1/m0/s1"
)
# Metal complexes: standard InChI removes the bonds to the metal, so a complex is written as its ligands and metal
# ions. Their largest fragment is a ligand, or the NH3 of any ammine complex, so they have no parent to compare.
CISPLATIN = "InChI=1S/2ClH.2H2N.Pt/h2*1H;2*1H2;/q;;2*-1;+4/p-2"  # N[Pt](N)(Cl)Cl
HEXAAMMINE_RUTHENIUM_CHLORIDE = (  # N[Ru](N)(N)(N)(N)N.[Cl-].[Cl-].[Cl-]
    "InChI=1S/3ClH.6H2N.Ru/h3*1H;6*1H2;/q;;;6*-1;+6/p-3"
)
OXALIPLATIN = (  # N[C@@H]1CCCC[C@H]1N.O=C1O[Pt]OC1=O
    "InChI=1S/C6H14N2.C2H2O4.Pt/c7-5-3-1-2-4-6(5)8;3-1(4)2(5)6;/h5-6H,1-4,7-8H2;(H,3,4)(H,5,6);/q;;+2/p-2"
    "/t5-,6-;;/m1../s1"
)
PT_DACH_DICHLORIDE = (  # N[C@@H]1CCCC[C@H]1N.Cl[Pt]Cl, the same diaminocyclohexane ligand as oxaliplatin
    "InChI=1S/C6H14N2.2ClH.Pt/c7-5-3-1-2-4-6(5)8;;;/h5-6H,1-4,7-8H2;2*1H;/q;;;+2/p-2/t5-,6-;;;/m1.../s1"
)
DACH = "InChI=1S/C6H14N2/c7-5-3-1-2-4-6(5)8/h5-6H,1-4,7-8H2/t5-,6-/m1/s1"  # N[C@@H]1CCCC[C@H]1N, the free ligand
# GtoPdb ligand 13351 (Ru360) and the PubChem record it links (CID 171663342): a bridging oxo against a water
RU360_GTOPDB = "InChI=1S/8H3N.2H2O.O.2Ru/h8*1H3;2*1H2;;;/q;;;;;;;;;;;+2;+3"
RU360_PUBCHEM = "InChI=1S/8H3N.3H2O.2Ru/h8*1H3;3*1H2;;/q;;;;;;;;;;;+1;+2/p+2"
YTTRIUM = "InChI=1S/Y"  # CHEBI:33331, which GtoPdb links to the peptide elamipretide
# the same formula, with a charge layer on one side only
SODIUM_ION = "InChI=1S/Na/q+1"  # [Na+]
SODIUM = "InChI=1S/Na"  # [Na]
ELAMIPRETIDE = (  # GtoPdb ligand 14276
    "InChI=1S/C32H49N9O5/c1-19-15-22(42)16-20(2)23(19)18-27(41-29(44)24(34)11-8-14-38-32(36)37)31(46)39-25(12-6-7-13"
    "-33)30(45)40-26(28(35)43)17-21-9-4-3-5-10-21/h3-5,9-10,15-16,24-27,42H,6-8,11-14,17-18,33-34H2,1-2H3,(H2,35,43)"
    "(H,39,46)(H,40,45)(H,41,44)(H4,36,37,38)/t24-,25+,26+,27+/m1/s1"
)
# GtoPdb ligand 2302 (ricasetron) and the PubChem record it links (CID 71785), which leaves the endo/exo centre 13
# undefined
RICASETRON_GTOPDB = (
    "InChI=1S/C19H27N3O/c1-19(2)12-22(17-7-5-4-6-16(17)19)18(23)20-13-10-14-8-9-15(11-13)21(14)3"
    "/h4-7,13-15H,8-12H2,1-3H3,(H,20,23)/t13-,14+,15-"
)
RICASETRON_PUBCHEM = (
    "InChI=1S/C19H27N3O/c1-19(2)12-22(17-7-5-4-6-16(17)19)18(23)20-13-10-14-8-9-15(11-13)21(14)3"
    "/h4-7,13-15H,8-12H2,1-3H3,(H,20,23)/t13?,14-,15+"
)
# GtoPdb ligand 14148 (AZ2158) and the ChEMBL record it links (CHEMBL1934424)
AZ2158_GTOPDB = (
    "InChI=1S/C25H38N4O4S/c1-16-8-14-21(15-9-16)34(31,32)28-23-22(24(30)26-18(3)25(4,5)6)17(2)27-29(23)19-10-12-"
    "20(33-7)13-11-19/h8-9,14-15,18-20,28H,10-13H2,1-7H3,(H,26,30)/t18-,19?,20-/m0/s1"
)
AZ2158_CHEMBL = (
    "InChI=1S/C25H38N4O4S/c1-16-8-14-21(15-9-16)34(31,32)28-23-22(24(30)26-18(3)25(4,5)6)17(2)27-29(23)19-10-12-"
    "20(33-7)13-11-19/h8-9,14-15,18-20,28H,10-13H2,1-7H3,(H,26,30)/t18-,19-,20+/m0/s1"
)
# calcium, like sodium, is a salt counterion rather than a metal complex, so these salts still have a parent
CALCIUM_ACETATE = "InChI=1S/2C2H4O2.Ca/c2*1-2(3)4;/h2*1H3,(H,3,4);/q;;+2/p-2"  # CC(=O)[O-].CC(=O)[O-].[Ca+2]
