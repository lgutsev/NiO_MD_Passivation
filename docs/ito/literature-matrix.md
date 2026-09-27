# ITO / SAM hole-selective contact literature matrix

- Search date: 2026-09-26
- Method: WebSearch + WebFetch against publisher pages, DOI resolvers, PubMed/PMC, arXiv/ChemRxiv. Rows added incrementally as verified.
- Verification legend: **V** = DOI and the cited details seen at source (full text or SI); **A** = abstract / publisher landing page only (details beyond the abstract are not confirmed); **U** = unverified (DOI or details not seen; what is missing is stated).
- Stack classes: `ITO/SAM` (SAM directly on ITO), `ITO/NiOx/SAM` (NiOx interlayer), `other` (FTO, ITO/Al2O3, glass, Au, etc.).

## Working log (in progress)

### [log] B3 verified 2026-09-26
- B3-1 Li D. et al., Nat. Commun. 15, 7605? (2024) "Co-adsorbed self-assembled monolayer enables high-performance perovskite and organic solar cells", DOI 10.1038/s41467-024-51760-5 (PMC11366757). Full text + SI read via Europe PMC / PMC / Springer ESM. V.
- B3-2 Park S.M. et al., Nature 624, 289-294 (2023), DOI 10.1038/s41586-023-06745-7. Accepted manuscript (Northwestern repository PDF) read; source of the MD recipe that Li et al. adapted. V (AAM).

### [log] A verified 2026-09-26
- A1 Al-Ashouri et al., EES 12, 3356 (2019), DOI 10.1039/C9EE02268F. Full text read from HZB repository copy (Unpaywall-listed). ITO/SAM, 0.5-1 mM EtOH (window 0.5-3 mM), spin or dip, 100 C 10 min, EtOH wash optional; RAIRS 1010 cm-1 P-O-ITO band + loss of P-OH ~950 cm-1 -> deprotonated anchor; WF ITO 4.6 eV, 2PACz 5.0 eV; DFT only for IR spectra and dipoles (2PACz +2 D, MeO-2PACz +0.2 D). V.
- A2 Al-Ashouri et al., Science 370, 1300 (2020), DOI 10.1126/science.abd4016. Author manuscript + SI read from HZB repository. ITO/SAM, 1 mM (~0.3 mg/ml) EtOH, 3000 rpm; ITO UV-O3 10-15 min before SAM; wash vs no-wash no device difference; dipoles ~1.7 D Me-4PACz; UPS WF shift in Fig S2 (values not transcribed). V.

### [log] A/B1/ITO-surface verified 2026-09-26
- Sindt C.A. et al. (Marder/Toney), ACS AMI 18, 29278 (2026), DOI 10.1021/acsami.6c03387 (PMC13220220). X-2PACz (H,F,Cl,Br,I,tBu) on UVO ITO and sapphire; 1 mM EtOH spin, 130 C 10 min, 3x EtOH rinse; NEXAFS carbazole tilt 61-65 deg; XRR coverage on sapphire 2PACz 2.31+/-0.38, I-2PACz 1.61+/-0.08 nm-2; ITO RMS 3.6 nm; no MD/DFT. V (main text via PMC; SI tables not read).
- Contreras H. et al. (Ginger/Armstrong), ACS AMI (2025), DOI 10.1021/acsami.5c12684 (DOE PAGES PDF read). I-2PACz 3 mM EtOH, spin (3000 rpm, 100 C 10 min) vs 12 h dip, O2 plasma 30 min vs HCl/FeCl3 etch; + 6dPA co-modifier; UPS WF table; O 1s hydroxyl 530.5-531.6 eV. V.
- ITO-surface refs harvested from Contreras reference list (DOIs seen in list; papers not opened unless noted): Brumbach 2007 Langmuir 10.1021/la701754u; Hotchkiss 2012 Acc Chem Res 10.1021/ar200119g; Paniagua 2008 JPCC 10.1021/jp710893k; Harrell 2018 JPCC 10.1021/acs.jpcc.7b10267; Chen 2022 ACS Nano 10.1021/acsnano.2c09115; Gliboff 2013 JPCC 10.1021/jp404033e.

