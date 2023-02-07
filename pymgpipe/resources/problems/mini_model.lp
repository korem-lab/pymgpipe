\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  communityBiomass - communityBiomass_reverse_5ad9d
Subject To
 acald[c]__TaxaA: - ACALD__TaxaA + ACALD__TaxaA_reverse_3bd4d
   + ACALDt__TaxaA - ACALDt__TaxaA_reverse_7366b + ALCD2x__TaxaA
   - ALCD2x__TaxaA_reverse_2678e = 0
 coa[c]__TaxaA: - ACALD__TaxaA + ACALD__TaxaA_reverse_3bd4d - AKGDH__TaxaA
   + AKGDH__TaxaA_reverse_48019 + 3.7478 biomass_ecoli__TaxaA
   - 3.7478 biomass_ecoli__TaxaA_reverse_ed191 + CS__TaxaA
   - CS__TaxaA_reverse_4b77c + MALS__TaxaA - MALS__TaxaA_reverse_c82e6
   - PDH__TaxaA + PDH__TaxaA_reverse_cbf23 - PFL__TaxaA
   + PFL__TaxaA_reverse_ff8bb + PTAr__TaxaA - PTAr__TaxaA_reverse_9e0a6
   - SUCOAS__TaxaA + SUCOAS__TaxaA_reverse_21627 = 0
 nad[c]__TaxaA: - ACALD__TaxaA + ACALD__TaxaA_reverse_3bd4d - AKGDH__TaxaA
   + AKGDH__TaxaA_reverse_48019 - ALCD2x__TaxaA
   + ALCD2x__TaxaA_reverse_2678e - 3.547 biomass_ecoli__TaxaA
   + 3.547 biomass_ecoli__TaxaA_reverse_ed191 - GAPD__TaxaA
   + GAPD__TaxaA_reverse_b06a4 - LDH_D__TaxaA + LDH_D__TaxaA_reverse_381f1
   - MDH__TaxaA + MDH__TaxaA_reverse_20572 - ME1__TaxaA
   + ME1__TaxaA_reverse_d8c3e + NADH16__TaxaA - NADH16__TaxaA_reverse_81d53
   - NADTRHD__TaxaA + NADTRHD__TaxaA_reverse_2a766 - PDH__TaxaA
   + PDH__TaxaA_reverse_cbf23 + THD2__TaxaA - THD2__TaxaA_reverse_f94fd
   = 0
 accoa[c]__TaxaA: ACALD__TaxaA - ACALD__TaxaA_reverse_3bd4d
   - 3.7478 biomass_ecoli__TaxaA
   + 3.7478 biomass_ecoli__TaxaA_reverse_ed191 - CS__TaxaA
   + CS__TaxaA_reverse_4b77c - MALS__TaxaA + MALS__TaxaA_reverse_c82e6
   + PDH__TaxaA - PDH__TaxaA_reverse_cbf23 + PFL__TaxaA
   - PFL__TaxaA_reverse_ff8bb - PTAr__TaxaA + PTAr__TaxaA_reverse_9e0a6
   = 0
 h[c]__TaxaA: ACALD__TaxaA - ACALD__TaxaA_reverse_3bd4d + ACt2r__TaxaA
   - ACt2r__TaxaA_reverse_5e5e9 + AKGt2r__TaxaA
   - AKGt2r__TaxaA_reverse_56a5c + ALCD2x__TaxaA
   - ALCD2x__TaxaA_reverse_2678e + ATPM__TaxaA - ATPM__TaxaA_reverse_d6dff
   + 3 ATPS4r__TaxaA - 3 ATPS4r__TaxaA_reverse_4a5d8
   + 59.81 biomass_ecoli__TaxaA - 59.81 biomass_ecoli__TaxaA_reverse_ed191
   + CS__TaxaA - CS__TaxaA_reverse_4b77c - 2 CYTBD__TaxaA
   + 2 CYTBD__TaxaA_reverse_983ef + D_LACt2__TaxaA
   - D_LACt2__TaxaA_reverse_a9e88 + ETOHt2r__TaxaA
   - ETOHt2r__TaxaA_reverse_8ca8d + FORt2__TaxaA
   - FORt2__TaxaA_reverse_7aa37 + 2 FUMt2_2__TaxaA
   - 2 FUMt2_2__TaxaA_reverse_ae205 + G6PDH2r__TaxaA
   - G6PDH2r__TaxaA_reverse_798d9 + GAPD__TaxaA - GAPD__TaxaA_reverse_b06a4
   + GLNS__TaxaA - GLNS__TaxaA_reverse_b05d2 + GLNabc__TaxaA
   - GLNabc__TaxaA_reverse_317f3 + GLUDy__TaxaA
   - GLUDy__TaxaA_reverse_fb35f - GLUSy__TaxaA + GLUSy__TaxaA_reverse_ec3ea
   + GLUt2r__TaxaA - GLUt2r__TaxaA_reverse_286c6 + LDH_D__TaxaA
   - LDH_D__TaxaA_reverse_381f1 + MALS__TaxaA - MALS__TaxaA_reverse_c82e6
   + 2 MALt2_2__TaxaA - 2 MALt2_2__TaxaA_reverse_87f6b + MDH__TaxaA
   - MDH__TaxaA_reverse_20572 - 4 NADH16__TaxaA
   + 4 NADH16__TaxaA_reverse_81d53 + PFK__TaxaA - PFK__TaxaA_reverse_ae474
   + PGL__TaxaA - PGL__TaxaA_reverse_23cba + PIt2r__TaxaA
   - PIt2r__TaxaA_reverse_e3f73 + PPC__TaxaA - PPC__TaxaA_reverse_39721
   + 2 PPS__TaxaA - 2 PPS__TaxaA_reverse_422f8 - PYK__TaxaA
   + PYK__TaxaA_reverse_dca3f + PYRt2__TaxaA - PYRt2__TaxaA_reverse_8f2e2
   + 2 SUCCt2_2__TaxaA - 2 SUCCt2_2__TaxaA_reverse_f20ba + SUCCt3__TaxaA
   - SUCCt3__TaxaA_reverse_f6546 + 2 THD2__TaxaA
   - 2 THD2__TaxaA_reverse_f94fd = 0
 nadh[c]__TaxaA: ACALD__TaxaA - ACALD__TaxaA_reverse_3bd4d + AKGDH__TaxaA
   - AKGDH__TaxaA_reverse_48019 + ALCD2x__TaxaA
   - ALCD2x__TaxaA_reverse_2678e + 3.547 biomass_ecoli__TaxaA
   - 3.547 biomass_ecoli__TaxaA_reverse_ed191 + GAPD__TaxaA
   - GAPD__TaxaA_reverse_b06a4 + LDH_D__TaxaA - LDH_D__TaxaA_reverse_381f1
   + MDH__TaxaA - MDH__TaxaA_reverse_20572 + ME1__TaxaA
   - ME1__TaxaA_reverse_d8c3e - NADH16__TaxaA + NADH16__TaxaA_reverse_81d53
   + NADTRHD__TaxaA - NADTRHD__TaxaA_reverse_2a766 + PDH__TaxaA
   - PDH__TaxaA_reverse_cbf23 - THD2__TaxaA + THD2__TaxaA_reverse_f94fd
   = 0
 acald[u]__TaxaA: - ACALDt__TaxaA + ACALDt__TaxaA_reverse_7366b
   + IEX_acald[u]__TaxaA - IEX_acald[u]__TaxaA_reverse_47467 = 0
 ac[c]__TaxaA: - ACKr__TaxaA + ACKr__TaxaA_reverse_35541 + ACt2r__TaxaA
   - ACt2r__TaxaA_reverse_5e5e9 = 0
 atp[c]__TaxaA: - ACKr__TaxaA + ACKr__TaxaA_reverse_35541 - ADK1__TaxaA
   + ADK1__TaxaA_reverse_af867 - ATPM__TaxaA + ATPM__TaxaA_reverse_d6dff
   + ATPS4r__TaxaA - ATPS4r__TaxaA_reverse_4a5d8
   - 59.81 biomass_ecoli__TaxaA + 59.81 biomass_ecoli__TaxaA_reverse_ed191
   - GLNS__TaxaA + GLNS__TaxaA_reverse_b05d2 - GLNabc__TaxaA
   + GLNabc__TaxaA_reverse_317f3 - PFK__TaxaA + PFK__TaxaA_reverse_ae474
   - PGK__TaxaA + PGK__TaxaA_reverse_3f06d - PPCK__TaxaA
   + PPCK__TaxaA_reverse_8d101 - PPS__TaxaA + PPS__TaxaA_reverse_422f8
   + PYK__TaxaA - PYK__TaxaA_reverse_dca3f - SUCOAS__TaxaA
   + SUCOAS__TaxaA_reverse_21627 = 0
 actp[c]__TaxaA: ACKr__TaxaA - ACKr__TaxaA_reverse_35541 + PTAr__TaxaA
   - PTAr__TaxaA_reverse_9e0a6 = 0
 adp[c]__TaxaA: ACKr__TaxaA - ACKr__TaxaA_reverse_35541 + 2 ADK1__TaxaA
   - 2 ADK1__TaxaA_reverse_af867 + ATPM__TaxaA - ATPM__TaxaA_reverse_d6dff
   - ATPS4r__TaxaA + ATPS4r__TaxaA_reverse_4a5d8
   + 59.81 biomass_ecoli__TaxaA - 59.81 biomass_ecoli__TaxaA_reverse_ed191
   + GLNS__TaxaA - GLNS__TaxaA_reverse_b05d2 + GLNabc__TaxaA
   - GLNabc__TaxaA_reverse_317f3 + PFK__TaxaA - PFK__TaxaA_reverse_ae474
   + PGK__TaxaA - PGK__TaxaA_reverse_3f06d + PPCK__TaxaA
   - PPCK__TaxaA_reverse_8d101 - PYK__TaxaA + PYK__TaxaA_reverse_dca3f
   + SUCOAS__TaxaA - SUCOAS__TaxaA_reverse_21627 = 0
 cit[c]__TaxaA: - ACONTa__TaxaA + ACONTa__TaxaA_reverse_dd81a + CS__TaxaA
   - CS__TaxaA_reverse_4b77c = 0
 acon_C[c]__TaxaA: ACONTa__TaxaA - ACONTa__TaxaA_reverse_dd81a
   - ACONTb__TaxaA + ACONTb__TaxaA_reverse_015ce = 0
 h2o[c]__TaxaA: ACONTa__TaxaA - ACONTa__TaxaA_reverse_dd81a - ACONTb__TaxaA
   + ACONTb__TaxaA_reverse_015ce - ATPM__TaxaA + ATPM__TaxaA_reverse_d6dff
   + ATPS4r__TaxaA - ATPS4r__TaxaA_reverse_4a5d8
   - 59.81 biomass_ecoli__TaxaA + 59.81 biomass_ecoli__TaxaA_reverse_ed191
   - CS__TaxaA + CS__TaxaA_reverse_4b77c + CYTBD__TaxaA
   - CYTBD__TaxaA_reverse_983ef + ENO__TaxaA - ENO__TaxaA_reverse_b3da1
   - FBP__TaxaA + FBP__TaxaA_reverse_8adae - FUM__TaxaA
   + FUM__TaxaA_reverse_9925b - GLNabc__TaxaA + GLNabc__TaxaA_reverse_317f3
   - GLUDy__TaxaA + GLUDy__TaxaA_reverse_fb35f - GLUN__TaxaA
   + GLUN__TaxaA_reverse_85542 + H2Ot__TaxaA - H2Ot__TaxaA_reverse_127eb
   - MALS__TaxaA + MALS__TaxaA_reverse_c82e6 - PGL__TaxaA
   + PGL__TaxaA_reverse_23cba - PPC__TaxaA + PPC__TaxaA_reverse_39721
   - PPS__TaxaA + PPS__TaxaA_reverse_422f8 = 0
 icit[c]__TaxaA: ACONTb__TaxaA - ACONTb__TaxaA_reverse_015ce
   - ICDHyr__TaxaA + ICDHyr__TaxaA_reverse_54fda - ICL__TaxaA
   + ICL__TaxaA_reverse_b0fb1 = 0
 ac[u]__TaxaA: - ACt2r__TaxaA + ACt2r__TaxaA_reverse_5e5e9
   + IEX_ac[u]__TaxaA - IEX_ac[u]__TaxaA_reverse_d85aa = 0
 h[u]__TaxaA: - ACt2r__TaxaA + ACt2r__TaxaA_reverse_5e5e9 - AKGt2r__TaxaA
   + AKGt2r__TaxaA_reverse_56a5c - 4 ATPS4r__TaxaA
   + 4 ATPS4r__TaxaA_reverse_4a5d8 + 2 CYTBD__TaxaA
   - 2 CYTBD__TaxaA_reverse_983ef - D_LACt2__TaxaA
   + D_LACt2__TaxaA_reverse_a9e88 - ETOHt2r__TaxaA
   + ETOHt2r__TaxaA_reverse_8ca8d + IEX_h[u]__TaxaA
   - IEX_h[u]__TaxaA_reverse_29989 - FORt2__TaxaA
   + FORt2__TaxaA_reverse_7aa37 - 2 FUMt2_2__TaxaA
   + 2 FUMt2_2__TaxaA_reverse_ae205 - GLUt2r__TaxaA
   + GLUt2r__TaxaA_reverse_286c6 - 2 MALt2_2__TaxaA
   + 2 MALt2_2__TaxaA_reverse_87f6b + 3 NADH16__TaxaA
   - 3 NADH16__TaxaA_reverse_81d53 - PIt2r__TaxaA
   + PIt2r__TaxaA_reverse_e3f73 - PYRt2__TaxaA + PYRt2__TaxaA_reverse_8f2e2
   - 2 SUCCt2_2__TaxaA + 2 SUCCt2_2__TaxaA_reverse_f20ba - SUCCt3__TaxaA
   + SUCCt3__TaxaA_reverse_f6546 - 2 THD2__TaxaA
   + 2 THD2__TaxaA_reverse_f94fd = 0
 amp[c]__TaxaA: - ADK1__TaxaA + ADK1__TaxaA_reverse_af867 + PPS__TaxaA
   - PPS__TaxaA_reverse_422f8 = 0
 akg[c]__TaxaA: - AKGDH__TaxaA + AKGDH__TaxaA_reverse_48019 + AKGt2r__TaxaA
   - AKGt2r__TaxaA_reverse_56a5c + 4.1182 biomass_ecoli__TaxaA
   - 4.1182 biomass_ecoli__TaxaA_reverse_ed191 + GLUDy__TaxaA
   - GLUDy__TaxaA_reverse_fb35f - GLUSy__TaxaA + GLUSy__TaxaA_reverse_ec3ea
   + ICDHyr__TaxaA - ICDHyr__TaxaA_reverse_54fda = 0
 co2[c]__TaxaA: AKGDH__TaxaA - AKGDH__TaxaA_reverse_48019 + CO2t__TaxaA
   - CO2t__TaxaA_reverse_5f188 + GND__TaxaA - GND__TaxaA_reverse_e35fe
   + ICDHyr__TaxaA - ICDHyr__TaxaA_reverse_54fda + ME1__TaxaA
   - ME1__TaxaA_reverse_d8c3e + ME2__TaxaA - ME2__TaxaA_reverse_b7ffb
   + PDH__TaxaA - PDH__TaxaA_reverse_cbf23 - PPC__TaxaA
   + PPC__TaxaA_reverse_39721 + PPCK__TaxaA - PPCK__TaxaA_reverse_8d101
   = 0
 succoa[c]__TaxaA: AKGDH__TaxaA - AKGDH__TaxaA_reverse_48019
   + SUCOAS__TaxaA - SUCOAS__TaxaA_reverse_21627 = 0
 akg[u]__TaxaA: - AKGt2r__TaxaA + AKGt2r__TaxaA_reverse_56a5c
   + IEX_akg[u]__TaxaA - IEX_akg[u]__TaxaA_reverse_21614 = 0
 etoh[c]__TaxaA: - ALCD2x__TaxaA + ALCD2x__TaxaA_reverse_2678e
   + ETOHt2r__TaxaA - ETOHt2r__TaxaA_reverse_8ca8d = 0
 pi[c]__TaxaA: ATPM__TaxaA - ATPM__TaxaA_reverse_d6dff - ATPS4r__TaxaA
   + ATPS4r__TaxaA_reverse_4a5d8 + 59.81 biomass_ecoli__TaxaA
   - 59.81 biomass_ecoli__TaxaA_reverse_ed191 + FBP__TaxaA
   - FBP__TaxaA_reverse_8adae - GAPD__TaxaA + GAPD__TaxaA_reverse_b06a4
   + GLNS__TaxaA - GLNS__TaxaA_reverse_b05d2 + GLNabc__TaxaA
   - GLNabc__TaxaA_reverse_317f3 + PIt2r__TaxaA
   - PIt2r__TaxaA_reverse_e3f73 + PPC__TaxaA - PPC__TaxaA_reverse_39721
   + PPS__TaxaA - PPS__TaxaA_reverse_422f8 - PTAr__TaxaA
   + PTAr__TaxaA_reverse_9e0a6 + SUCOAS__TaxaA
   - SUCOAS__TaxaA_reverse_21627 = 0
 3pg[c]__TaxaA: - 1.496 biomass_ecoli__TaxaA
   + 1.496 biomass_ecoli__TaxaA_reverse_ed191 - PGK__TaxaA
   + PGK__TaxaA_reverse_3f06d + PGM__TaxaA - PGM__TaxaA_reverse_b2abe = 0
 e4p[c]__TaxaA: - 0.361 biomass_ecoli__TaxaA
   + 0.361 biomass_ecoli__TaxaA_reverse_ed191 + TALA__TaxaA
   - TALA__TaxaA_reverse_0bcac - TKT2__TaxaA + TKT2__TaxaA_reverse_779d3
   = 0
 f6p[c]__TaxaA: - 0.0709 biomass_ecoli__TaxaA
   + 0.0709 biomass_ecoli__TaxaA_reverse_ed191 + FBP__TaxaA
   - FBP__TaxaA_reverse_8adae + FRUpts2__TaxaA
   - FRUpts2__TaxaA_reverse_30526 - PFK__TaxaA + PFK__TaxaA_reverse_ae474
   + PGI__TaxaA - PGI__TaxaA_reverse_c02d7 + TALA__TaxaA
   - TALA__TaxaA_reverse_0bcac + TKT2__TaxaA - TKT2__TaxaA_reverse_779d3
   = 0
 g3p[c]__TaxaA: - 0.129 biomass_ecoli__TaxaA
   + 0.129 biomass_ecoli__TaxaA_reverse_ed191 + FBA__TaxaA
   - FBA__TaxaA_reverse_7aab1 - GAPD__TaxaA + GAPD__TaxaA_reverse_b06a4
   - TALA__TaxaA + TALA__TaxaA_reverse_0bcac + TKT1__TaxaA
   - TKT1__TaxaA_reverse_3bb31 + TKT2__TaxaA - TKT2__TaxaA_reverse_779d3
   + TPI__TaxaA - TPI__TaxaA_reverse_9fc16 = 0
 g6p[c]__TaxaA: - 0.205 biomass_ecoli__TaxaA
   + 0.205 biomass_ecoli__TaxaA_reverse_ed191 - G6PDH2r__TaxaA
   + G6PDH2r__TaxaA_reverse_798d9 + GLCpts__TaxaA
   - GLCpts__TaxaA_reverse_d6498 - PGI__TaxaA + PGI__TaxaA_reverse_c02d7
   = 0
 gln__L[c]__TaxaA: - 0.2557 biomass_ecoli__TaxaA
   + 0.2557 biomass_ecoli__TaxaA_reverse_ed191 + GLNS__TaxaA
   - GLNS__TaxaA_reverse_b05d2 + GLNabc__TaxaA
   - GLNabc__TaxaA_reverse_317f3 - GLUN__TaxaA + GLUN__TaxaA_reverse_85542
   - GLUSy__TaxaA + GLUSy__TaxaA_reverse_ec3ea = 0
 glu__L[c]__TaxaA: - 4.9414 biomass_ecoli__TaxaA
   + 4.9414 biomass_ecoli__TaxaA_reverse_ed191 - GLNS__TaxaA
   + GLNS__TaxaA_reverse_b05d2 - GLUDy__TaxaA + GLUDy__TaxaA_reverse_fb35f
   + GLUN__TaxaA - GLUN__TaxaA_reverse_85542 + 2 GLUSy__TaxaA
   - 2 GLUSy__TaxaA_reverse_ec3ea + GLUt2r__TaxaA
   - GLUt2r__TaxaA_reverse_286c6 = 0
 nadph[c]__TaxaA: - 13.0279 biomass_ecoli__TaxaA
   + 13.0279 biomass_ecoli__TaxaA_reverse_ed191 + G6PDH2r__TaxaA
   - G6PDH2r__TaxaA_reverse_798d9 + GLUDy__TaxaA
   - GLUDy__TaxaA_reverse_fb35f - GLUSy__TaxaA + GLUSy__TaxaA_reverse_ec3ea
   + GND__TaxaA - GND__TaxaA_reverse_e35fe + ICDHyr__TaxaA
   - ICDHyr__TaxaA_reverse_54fda + ME2__TaxaA - ME2__TaxaA_reverse_b7ffb
   - NADTRHD__TaxaA + NADTRHD__TaxaA_reverse_2a766 + THD2__TaxaA
   - THD2__TaxaA_reverse_f94fd = 0
 oaa[c]__TaxaA: - 1.7867 biomass_ecoli__TaxaA
   + 1.7867 biomass_ecoli__TaxaA_reverse_ed191 - CS__TaxaA
   + CS__TaxaA_reverse_4b77c + MDH__TaxaA - MDH__TaxaA_reverse_20572
   + PPC__TaxaA - PPC__TaxaA_reverse_39721 - PPCK__TaxaA
   + PPCK__TaxaA_reverse_8d101 = 0
 pep[c]__TaxaA: - 0.5191 biomass_ecoli__TaxaA
   + 0.5191 biomass_ecoli__TaxaA_reverse_ed191 + ENO__TaxaA
   - ENO__TaxaA_reverse_b3da1 - FRUpts2__TaxaA
   + FRUpts2__TaxaA_reverse_30526 - GLCpts__TaxaA
   + GLCpts__TaxaA_reverse_d6498 - PPC__TaxaA + PPC__TaxaA_reverse_39721
   + PPCK__TaxaA - PPCK__TaxaA_reverse_8d101 + PPS__TaxaA
   - PPS__TaxaA_reverse_422f8 - PYK__TaxaA + PYK__TaxaA_reverse_dca3f = 0
 pyr[c]__TaxaA: - 2.8328 biomass_ecoli__TaxaA
   + 2.8328 biomass_ecoli__TaxaA_reverse_ed191 + FRUpts2__TaxaA
   - FRUpts2__TaxaA_reverse_30526 + GLCpts__TaxaA
   - GLCpts__TaxaA_reverse_d6498 + LDH_D__TaxaA
   - LDH_D__TaxaA_reverse_381f1 + ME1__TaxaA - ME1__TaxaA_reverse_d8c3e
   + ME2__TaxaA - ME2__TaxaA_reverse_b7ffb - PDH__TaxaA
   + PDH__TaxaA_reverse_cbf23 - PFL__TaxaA + PFL__TaxaA_reverse_ff8bb
   - PPS__TaxaA + PPS__TaxaA_reverse_422f8 + PYK__TaxaA
   - PYK__TaxaA_reverse_dca3f + PYRt2__TaxaA - PYRt2__TaxaA_reverse_8f2e2
   = 0
 r5p[c]__TaxaA: - 0.8977 biomass_ecoli__TaxaA
   + 0.8977 biomass_ecoli__TaxaA_reverse_ed191 - RPI__TaxaA
   + RPI__TaxaA_reverse_b8683 - TKT1__TaxaA + TKT1__TaxaA_reverse_3bb31
   = 0
 nadp[c]__TaxaA: 13.0279 biomass_ecoli__TaxaA
   - 13.0279 biomass_ecoli__TaxaA_reverse_ed191 - G6PDH2r__TaxaA
   + G6PDH2r__TaxaA_reverse_798d9 - GLUDy__TaxaA
   + GLUDy__TaxaA_reverse_fb35f + GLUSy__TaxaA - GLUSy__TaxaA_reverse_ec3ea
   - GND__TaxaA + GND__TaxaA_reverse_e35fe - ICDHyr__TaxaA
   + ICDHyr__TaxaA_reverse_54fda - ME2__TaxaA + ME2__TaxaA_reverse_b7ffb
   + NADTRHD__TaxaA - NADTRHD__TaxaA_reverse_2a766 - THD2__TaxaA
   + THD2__TaxaA_reverse_f94fd = 0
 biomass[c]__TaxaA: biomass_ecoli__TaxaA
   - biomass_ecoli__TaxaA_reverse_ed191 - 0.1 communityBiomass
   + 0.1 communityBiomass_reverse_5ad9d = 0
 co2[u]__TaxaA: - CO2t__TaxaA + CO2t__TaxaA_reverse_5f188
   + IEX_co2[u]__TaxaA - IEX_co2[u]__TaxaA_reverse_3ae77 = 0
 o2[c]__TaxaA: - 0.5 CYTBD__TaxaA + 0.5 CYTBD__TaxaA_reverse_983ef
   + O2t__TaxaA - O2t__TaxaA_reverse_b8a83 = 0
 q8h2[c]__TaxaA: - CYTBD__TaxaA + CYTBD__TaxaA_reverse_983ef - FRD7__TaxaA
   + FRD7__TaxaA_reverse_80d05 + NADH16__TaxaA
   - NADH16__TaxaA_reverse_81d53 + SUCDi__TaxaA
   - SUCDi__TaxaA_reverse_6fb91 = 0
 q8[c]__TaxaA: CYTBD__TaxaA - CYTBD__TaxaA_reverse_983ef + FRD7__TaxaA
   - FRD7__TaxaA_reverse_80d05 - NADH16__TaxaA
   + NADH16__TaxaA_reverse_81d53 - SUCDi__TaxaA
   + SUCDi__TaxaA_reverse_6fb91 = 0
 lac__D[u]__TaxaA: - D_LACt2__TaxaA + D_LACt2__TaxaA_reverse_a9e88
   + IEX_lac__D[u]__TaxaA - IEX_lac__D[u]__TaxaA_reverse_d234a = 0
 lac__D[c]__TaxaA: D_LACt2__TaxaA - D_LACt2__TaxaA_reverse_a9e88
   - LDH_D__TaxaA + LDH_D__TaxaA_reverse_381f1 = 0
 2pg[c]__TaxaA: - ENO__TaxaA + ENO__TaxaA_reverse_b3da1 - PGM__TaxaA
   + PGM__TaxaA_reverse_b2abe = 0
 etoh[u]__TaxaA: - ETOHt2r__TaxaA + ETOHt2r__TaxaA_reverse_8ca8d
   + EX[u]toh[u]__TaxaA - EX[u]toh[u]__TaxaA_reverse_c58a3 = 0
 for[u]__TaxaA: IEX_for[u]__TaxaA - IEX_for[u]__TaxaA_reverse_f0688
   - FORt2__TaxaA + FORt2__TaxaA_reverse_7aa37 + FORti__TaxaA
   - FORti__TaxaA_reverse_2b70c = 0
 fru[u]__TaxaA: IEX_fru[u]__TaxaA - IEX_fru[u]__TaxaA_reverse_0ac79
   - FRUpts2__TaxaA + FRUpts2__TaxaA_reverse_30526 = 0
 fum[u]__TaxaA: IEX_fum[u]__TaxaA - IEX_fum[u]__TaxaA_reverse_b81bc
   - FUMt2_2__TaxaA + FUMt2_2__TaxaA_reverse_ae205 = 0
 glc__D[u]__TaxaA: IEX_glc__D[u]__TaxaA
   - IEX_glc__D[u]__TaxaA_reverse_e84a1 - GLCpts__TaxaA
   + GLCpts__TaxaA_reverse_d6498 = 0
 gln__L[u]__TaxaA: IEX_gln__L[u]__TaxaA
   - IEX_gln__L[u]__TaxaA_reverse_afdbf - GLNabc__TaxaA
   + GLNabc__TaxaA_reverse_317f3 = 0
 glu__L[u]__TaxaA: IEX_glu__L[u]__TaxaA
   - IEX_glu__L[u]__TaxaA_reverse_2873c - GLUt2r__TaxaA
   + GLUt2r__TaxaA_reverse_286c6 = 0
 h2o[u]__TaxaA: IEX_h2o[u]__TaxaA - IEX_h2o[u]__TaxaA_reverse_dd098
   - H2Ot__TaxaA + H2Ot__TaxaA_reverse_127eb = 0
 mal__L[u]__TaxaA: IEX_mal__L[u]__TaxaA
   - IEX_mal__L[u]__TaxaA_reverse_f589a - MALt2_2__TaxaA
   + MALt2_2__TaxaA_reverse_87f6b = 0
 nh4[u]__TaxaA: IEX_nh4[u]__TaxaA - IEX_nh4[u]__TaxaA_reverse_d1e33
   - NH4t__TaxaA + NH4t__TaxaA_reverse_7a63f = 0
 o2[u]__TaxaA: IEX_o2[u]__TaxaA - IEX_o2[u]__TaxaA_reverse_cbfe4
   - O2t__TaxaA + O2t__TaxaA_reverse_b8a83 = 0
 pi[u]__TaxaA: IEX_pi[u]__TaxaA - IEX_pi[u]__TaxaA_reverse_621ea
   - PIt2r__TaxaA + PIt2r__TaxaA_reverse_e3f73 = 0
 pyr[u]__TaxaA: IEX_pyr[u]__TaxaA - IEX_pyr[u]__TaxaA_reverse_07e94
   - PYRt2__TaxaA + PYRt2__TaxaA_reverse_8f2e2 = 0
 succ[u]__TaxaA: IEX_succ[u]__TaxaA - IEX_succ[u]__TaxaA_reverse_75565
   - SUCCt2_2__TaxaA + SUCCt2_2__TaxaA_reverse_f20ba + SUCCt3__TaxaA
   - SUCCt3__TaxaA_reverse_f6546 = 0
 fdp[c]__TaxaA: - FBA__TaxaA + FBA__TaxaA_reverse_7aab1 - FBP__TaxaA
   + FBP__TaxaA_reverse_8adae + PFK__TaxaA - PFK__TaxaA_reverse_ae474 = 0
 dhap[c]__TaxaA: FBA__TaxaA - FBA__TaxaA_reverse_7aab1 - TPI__TaxaA
   + TPI__TaxaA_reverse_9fc16 = 0
 for[c]__TaxaA: FORt2__TaxaA - FORt2__TaxaA_reverse_7aa37 - FORti__TaxaA
   + FORti__TaxaA_reverse_2b70c + PFL__TaxaA - PFL__TaxaA_reverse_ff8bb
   = 0
 fum[c]__TaxaA: - FRD7__TaxaA + FRD7__TaxaA_reverse_80d05 - FUM__TaxaA
   + FUM__TaxaA_reverse_9925b + FUMt2_2__TaxaA
   - FUMt2_2__TaxaA_reverse_ae205 + SUCDi__TaxaA
   - SUCDi__TaxaA_reverse_6fb91 = 0
 succ[c]__TaxaA: FRD7__TaxaA - FRD7__TaxaA_reverse_80d05 + ICL__TaxaA
   - ICL__TaxaA_reverse_b0fb1 + SUCCt2_2__TaxaA
   - SUCCt2_2__TaxaA_reverse_f20ba - SUCCt3__TaxaA
   + SUCCt3__TaxaA_reverse_f6546 - SUCDi__TaxaA
   + SUCDi__TaxaA_reverse_6fb91 - SUCOAS__TaxaA
   + SUCOAS__TaxaA_reverse_21627 = 0
 mal__L[c]__TaxaA: FUM__TaxaA - FUM__TaxaA_reverse_9925b + MALS__TaxaA
   - MALS__TaxaA_reverse_c82e6 + MALt2_2__TaxaA
   - MALt2_2__TaxaA_reverse_87f6b - MDH__TaxaA + MDH__TaxaA_reverse_20572
   - ME1__TaxaA + ME1__TaxaA_reverse_d8c3e - ME2__TaxaA
   + ME2__TaxaA_reverse_b7ffb = 0
 6pgl[c]__TaxaA: G6PDH2r__TaxaA - G6PDH2r__TaxaA_reverse_798d9 - PGL__TaxaA
   + PGL__TaxaA_reverse_23cba = 0
 13dpg[c]__TaxaA: GAPD__TaxaA - GAPD__TaxaA_reverse_b06a4 + PGK__TaxaA
   - PGK__TaxaA_reverse_3f06d = 0
 nh4[c]__TaxaA: - GLNS__TaxaA + GLNS__TaxaA_reverse_b05d2 + GLUDy__TaxaA
   - GLUDy__TaxaA_reverse_fb35f + GLUN__TaxaA - GLUN__TaxaA_reverse_85542
   + NH4t__TaxaA - NH4t__TaxaA_reverse_7a63f = 0
 6pgc[c]__TaxaA: - GND__TaxaA + GND__TaxaA_reverse_e35fe + PGL__TaxaA
   - PGL__TaxaA_reverse_23cba = 0
 ru5p__D[c]__TaxaA: GND__TaxaA - GND__TaxaA_reverse_e35fe - RPE__TaxaA
   + RPE__TaxaA_reverse_50231 + RPI__TaxaA - RPI__TaxaA_reverse_b8683 = 0
 glx[c]__TaxaA: ICL__TaxaA - ICL__TaxaA_reverse_b0fb1 - MALS__TaxaA
   + MALS__TaxaA_reverse_c82e6 = 0
 xu5p__D[c]__TaxaA: RPE__TaxaA - RPE__TaxaA_reverse_50231 - TKT1__TaxaA
   + TKT1__TaxaA_reverse_3bb31 - TKT2__TaxaA + TKT2__TaxaA_reverse_779d3
   = 0
 s7p[c]__TaxaA: - TALA__TaxaA + TALA__TaxaA_reverse_0bcac + TKT1__TaxaA
   - TKT1__TaxaA_reverse_3bb31 = 0
 ac[u]: - IEX_ac[u]__TaxaA + IEX_ac[u]__TaxaA_reverse_d85aa - UFEt_ac
   + UFEt_ac_reverse_0351f + DUt_ac - DUt_ac_reverse_d6d71
   - IEX_ac[u]__TaxaB + IEX_ac[u]__TaxaB_reverse_b168f - IEX_ac[u]__TaxaC
   + IEX_ac[u]__TaxaC_reverse_82616 - IEX_ac[u]__TaxaD
   + IEX_ac[u]__TaxaD_reverse_63ca5 = 0
 ac[fe]: - EX_ac[fe] + EX_ac[fe]_reverse_15dd0 + UFEt_ac
   - UFEt_ac_reverse_0351f = 0
 ac[d]: - Diet_EX_ac[d] + Diet_EX_ac[d]_reverse_c5b1c - DUt_ac
   + DUt_ac_reverse_d6d71 = 0
 acald[u]: - IEX_acald[u]__TaxaA + IEX_acald[u]__TaxaA_reverse_47467
   - UFEt_acald + UFEt_acald_reverse_23948 + DUt_acald
   - DUt_acald_reverse_88448 - IEX_acald[u]__TaxaB
   + IEX_acald[u]__TaxaB_reverse_2ec49 - IEX_acald[u]__TaxaC
   + IEX_acald[u]__TaxaC_reverse_79433 - IEX_acald[u]__TaxaD
   + IEX_acald[u]__TaxaD_reverse_1bc86 = 0
 acald[fe]: - EX_acald[fe] + EX_acald[fe]_reverse_745af + UFEt_acald
   - UFEt_acald_reverse_23948 = 0
 acald[d]: - Diet_EX_acald[d] + Diet_EX_acald[d]_reverse_dcfe5 - DUt_acald
   + DUt_acald_reverse_88448 = 0
 akg[u]: - IEX_akg[u]__TaxaA + IEX_akg[u]__TaxaA_reverse_21614 - UFEt_akg
   + UFEt_akg_reverse_a80ce + DUt_akg - DUt_akg_reverse_eb159
   - IEX_akg[u]__TaxaB + IEX_akg[u]__TaxaB_reverse_427b9
   - IEX_akg[u]__TaxaC + IEX_akg[u]__TaxaC_reverse_2d021
   - IEX_akg[u]__TaxaD + IEX_akg[u]__TaxaD_reverse_a37cb = 0
 akg[fe]: - EX_akg[fe] + EX_akg[fe]_reverse_3e084 + UFEt_akg
   - UFEt_akg_reverse_a80ce = 0
 akg[d]: - Diet_EX_akg[d] + Diet_EX_akg[d]_reverse_e3523 - DUt_akg
   + DUt_akg_reverse_eb159 = 0
 co2[u]: - IEX_co2[u]__TaxaA + IEX_co2[u]__TaxaA_reverse_3ae77 - UFEt_co2
   + UFEt_co2_reverse_de6bf + DUt_co2 - DUt_co2_reverse_192a6
   - IEX_co2[u]__TaxaB + IEX_co2[u]__TaxaB_reverse_60e91
   - IEX_co2[u]__TaxaC + IEX_co2[u]__TaxaC_reverse_f8836
   - IEX_co2[u]__TaxaD + IEX_co2[u]__TaxaD_reverse_ce75e = 0
 co2[fe]: - EX_co2[fe] + EX_co2[fe]_reverse_9f050 + UFEt_co2
   - UFEt_co2_reverse_de6bf = 0
 co2[d]: - Diet_EX_co2[d] + Diet_EX_co2[d]_reverse_ec81a - DUt_co2
   + DUt_co2_reverse_192a6 = 0
 etoh[u]: - EX[u]toh[u]__TaxaA + EX[u]toh[u]__TaxaA_reverse_c58a3
   - UFEt_etoh + UFEt_etoh_reverse_bbee7 + DUt_etoh
   - DUt_etoh_reverse_c3f32 - EX[u]toh[u]__TaxaB
   + EX[u]toh[u]__TaxaB_reverse_0c398 - EX[u]toh[u]__TaxaC
   + EX[u]toh[u]__TaxaC_reverse_65089 - EX[u]toh[u]__TaxaD
   + EX[u]toh[u]__TaxaD_reverse_1d69d = 0
 etoh[fe]: - EX_etoh[fe] + EX_etoh[fe]_reverse_406f6 + UFEt_etoh
   - UFEt_etoh_reverse_bbee7 = 0
 etoh[d]: - Diet_EX_etoh[d] + Diet_EX_etoh[d]_reverse_f4854 - DUt_etoh
   + DUt_etoh_reverse_c3f32 = 0
 for[u]: - IEX_for[u]__TaxaA + IEX_for[u]__TaxaA_reverse_f0688 - UFEt_for
   + UFEt_for_reverse_7662b + DUt_for - DUt_for_reverse_f8c5e
   - IEX_for[u]__TaxaB + IEX_for[u]__TaxaB_reverse_10439
   - IEX_for[u]__TaxaC + IEX_for[u]__TaxaC_reverse_89537
   - IEX_for[u]__TaxaD + IEX_for[u]__TaxaD_reverse_0a866 = 0
 for[fe]: - EX_for[fe] + EX_for[fe]_reverse_5e89a + UFEt_for
   - UFEt_for_reverse_7662b = 0
 for[d]: - Diet_EX_for[d] + Diet_EX_for[d]_reverse_4b919 - DUt_for
   + DUt_for_reverse_f8c5e = 0
 fru[u]: - IEX_fru[u]__TaxaA + IEX_fru[u]__TaxaA_reverse_0ac79 - UFEt_fru
   + UFEt_fru_reverse_1d253 + DUt_fru - DUt_fru_reverse_b898a
   - IEX_fru[u]__TaxaB + IEX_fru[u]__TaxaB_reverse_3b7e3
   - IEX_fru[u]__TaxaC + IEX_fru[u]__TaxaC_reverse_d1617
   - IEX_fru[u]__TaxaD + IEX_fru[u]__TaxaD_reverse_8c9bf = 0
 fru[fe]: - EX_fru[fe] + EX_fru[fe]_reverse_9e214 + UFEt_fru
   - UFEt_fru_reverse_1d253 = 0
 fru[d]: - Diet_EX_fru[d] + Diet_EX_fru[d]_reverse_079dc - DUt_fru
   + DUt_fru_reverse_b898a = 0
 fum[u]: - IEX_fum[u]__TaxaA + IEX_fum[u]__TaxaA_reverse_b81bc - UFEt_fum
   + UFEt_fum_reverse_9d01c + DUt_fum - DUt_fum_reverse_89aaf
   - IEX_fum[u]__TaxaB + IEX_fum[u]__TaxaB_reverse_14c2d
   - IEX_fum[u]__TaxaC + IEX_fum[u]__TaxaC_reverse_66a1c
   - IEX_fum[u]__TaxaD + IEX_fum[u]__TaxaD_reverse_3415d = 0
 fum[fe]: - EX_fum[fe] + EX_fum[fe]_reverse_7136c + UFEt_fum
   - UFEt_fum_reverse_9d01c = 0
 fum[d]: - Diet_EX_fum[d] + Diet_EX_fum[d]_reverse_818ef - DUt_fum
   + DUt_fum_reverse_89aaf = 0
 glc__D[u]: - IEX_glc__D[u]__TaxaA + IEX_glc__D[u]__TaxaA_reverse_e84a1
   - UFEt_glc__D + UFEt_glc__D_reverse_691ca + DUt_glc__D
   - DUt_glc__D_reverse_c951b - IEX_glc__D[u]__TaxaB
   + IEX_glc__D[u]__TaxaB_reverse_37fbc - IEX_glc__D[u]__TaxaC
   + IEX_glc__D[u]__TaxaC_reverse_4dbc9 - IEX_glc__D[u]__TaxaD
   + IEX_glc__D[u]__TaxaD_reverse_405dc = 0
 glc__D[fe]: - EX_glc__D[fe] + EX_glc__D[fe]_reverse_48473 + UFEt_glc__D
   - UFEt_glc__D_reverse_691ca = 0
 glc__D[d]: - Diet_EX_glc__D[d] + Diet_EX_glc__D[d]_reverse_a4b1e
   - DUt_glc__D + DUt_glc__D_reverse_c951b = 0
 gln__L[u]: - IEX_gln__L[u]__TaxaA + IEX_gln__L[u]__TaxaA_reverse_afdbf
   - UFEt_gln__L + UFEt_gln__L_reverse_76d9b + DUt_gln__L
   - DUt_gln__L_reverse_374f9 - IEX_gln__L[u]__TaxaB
   + IEX_gln__L[u]__TaxaB_reverse_eb617 - IEX_gln__L[u]__TaxaC
   + IEX_gln__L[u]__TaxaC_reverse_d4f65 - IEX_gln__L[u]__TaxaD
   + IEX_gln__L[u]__TaxaD_reverse_23702 = 0
 gln__L[fe]: - EX_gln__L[fe] + EX_gln__L[fe]_reverse_f5e4f + UFEt_gln__L
   - UFEt_gln__L_reverse_76d9b = 0
 gln__L[d]: - Diet_EX_gln__L[d] + Diet_EX_gln__L[d]_reverse_2ef88
   - DUt_gln__L + DUt_gln__L_reverse_374f9 = 0
 glu__L[u]: - IEX_glu__L[u]__TaxaA + IEX_glu__L[u]__TaxaA_reverse_2873c
   - UFEt_glu__L + UFEt_glu__L_reverse_e0c79 + DUt_glu__L
   - DUt_glu__L_reverse_331f7 - IEX_glu__L[u]__TaxaB
   + IEX_glu__L[u]__TaxaB_reverse_afc30 - IEX_glu__L[u]__TaxaC
   + IEX_glu__L[u]__TaxaC_reverse_b053a - IEX_glu__L[u]__TaxaD
   + IEX_glu__L[u]__TaxaD_reverse_db289 = 0
 glu__L[fe]: - EX_glu__L[fe] + EX_glu__L[fe]_reverse_24640 + UFEt_glu__L
   - UFEt_glu__L_reverse_e0c79 = 0
 glu__L[d]: - Diet_EX_glu__L[d] + Diet_EX_glu__L[d]_reverse_6886b
   - DUt_glu__L + DUt_glu__L_reverse_331f7 = 0
 h[u]: - IEX_h[u]__TaxaA + IEX_h[u]__TaxaA_reverse_29989 - UFEt_h
   + UFEt_h_reverse_0e86f + DUt_h - DUt_h_reverse_e9a4a - IEX_h[u]__TaxaB
   + IEX_h[u]__TaxaB_reverse_a4431 - IEX_h[u]__TaxaC
   + IEX_h[u]__TaxaC_reverse_b39fd - IEX_h[u]__TaxaD
   + IEX_h[u]__TaxaD_reverse_3e6c7 = 0
 h[fe]: - EX_h[fe] + EX_h[fe]_reverse_f5bf8 + UFEt_h - UFEt_h_reverse_0e86f
   = 0
 h[d]: - Diet_EX_h[d] + Diet_EX_h[d]_reverse_04d8d - DUt_h
   + DUt_h_reverse_e9a4a = 0
 h2o[u]: - IEX_h2o[u]__TaxaA + IEX_h2o[u]__TaxaA_reverse_dd098 - UFEt_h2o
   + UFEt_h2o_reverse_9fc26 + DUt_h2o - DUt_h2o_reverse_067ac
   - IEX_h2o[u]__TaxaB + IEX_h2o[u]__TaxaB_reverse_32bd0
   - IEX_h2o[u]__TaxaC + IEX_h2o[u]__TaxaC_reverse_f3c2c
   - IEX_h2o[u]__TaxaD + IEX_h2o[u]__TaxaD_reverse_c1c75 = 0
 h2o[fe]: - EX_h2o[fe] + EX_h2o[fe]_reverse_e4723 + UFEt_h2o
   - UFEt_h2o_reverse_9fc26 = 0
 h2o[d]: - Diet_EX_h2o[d] + Diet_EX_h2o[d]_reverse_1812e - DUt_h2o
   + DUt_h2o_reverse_067ac = 0
 lac__D[u]: - IEX_lac__D[u]__TaxaA + IEX_lac__D[u]__TaxaA_reverse_d234a
   - UFEt_lac__D + UFEt_lac__D_reverse_7be71 + DUt_lac__D
   - DUt_lac__D_reverse_cb570 - IEX_lac__D[u]__TaxaB
   + IEX_lac__D[u]__TaxaB_reverse_03348 - IEX_lac__D[u]__TaxaC
   + IEX_lac__D[u]__TaxaC_reverse_baf61 - IEX_lac__D[u]__TaxaD
   + IEX_lac__D[u]__TaxaD_reverse_115c3 = 0
 lac__D[fe]: - EX_lac__D[fe] + EX_lac__D[fe]_reverse_57587 + UFEt_lac__D
   - UFEt_lac__D_reverse_7be71 = 0
 lac__D[d]: - Diet_EX_lac__D[d] + Diet_EX_lac__D[d]_reverse_07348
   - DUt_lac__D + DUt_lac__D_reverse_cb570 = 0
 mal__L[u]: - IEX_mal__L[u]__TaxaA + IEX_mal__L[u]__TaxaA_reverse_f589a
   - UFEt_mal__L + UFEt_mal__L_reverse_89431 + DUt_mal__L
   - DUt_mal__L_reverse_4cd47 - IEX_mal__L[u]__TaxaB
   + IEX_mal__L[u]__TaxaB_reverse_e69bf - IEX_mal__L[u]__TaxaC
   + IEX_mal__L[u]__TaxaC_reverse_2b6d1 - IEX_mal__L[u]__TaxaD
   + IEX_mal__L[u]__TaxaD_reverse_36f62 = 0
 mal__L[fe]: - EX_mal__L[fe] + EX_mal__L[fe]_reverse_82659 + UFEt_mal__L
   - UFEt_mal__L_reverse_89431 = 0
 mal__L[d]: - Diet_EX_mal__L[d] + Diet_EX_mal__L[d]_reverse_44e05
   - DUt_mal__L + DUt_mal__L_reverse_4cd47 = 0
 nh4[u]: - IEX_nh4[u]__TaxaA + IEX_nh4[u]__TaxaA_reverse_d1e33 - UFEt_nh4
   + UFEt_nh4_reverse_f9f09 + DUt_nh4 - DUt_nh4_reverse_1b2f2
   - IEX_nh4[u]__TaxaB + IEX_nh4[u]__TaxaB_reverse_e819b
   - IEX_nh4[u]__TaxaC + IEX_nh4[u]__TaxaC_reverse_486f8
   - IEX_nh4[u]__TaxaD + IEX_nh4[u]__TaxaD_reverse_a85a3 = 0
 nh4[fe]: - EX_nh4[fe] + EX_nh4[fe]_reverse_02041 + UFEt_nh4
   - UFEt_nh4_reverse_f9f09 = 0
 nh4[d]: - Diet_EX_nh4[d] + Diet_EX_nh4[d]_reverse_4a802 - DUt_nh4
   + DUt_nh4_reverse_1b2f2 = 0
 o2[u]: - IEX_o2[u]__TaxaA + IEX_o2[u]__TaxaA_reverse_cbfe4 - UFEt_o2
   + UFEt_o2_reverse_d47cc + DUt_o2 - DUt_o2_reverse_5b3e1
   - IEX_o2[u]__TaxaB + IEX_o2[u]__TaxaB_reverse_6c7da - IEX_o2[u]__TaxaC
   + IEX_o2[u]__TaxaC_reverse_9d59a - IEX_o2[u]__TaxaD
   + IEX_o2[u]__TaxaD_reverse_a1b64 = 0
 o2[fe]: - EX_o2[fe] + EX_o2[fe]_reverse_4a084 + UFEt_o2
   - UFEt_o2_reverse_d47cc = 0
 o2[d]: - Diet_EX_o2[d] + Diet_EX_o2[d]_reverse_6f9b0 - DUt_o2
   + DUt_o2_reverse_5b3e1 = 0
 pi[u]: - IEX_pi[u]__TaxaA + IEX_pi[u]__TaxaA_reverse_621ea - UFEt_pi
   + UFEt_pi_reverse_d3357 + DUt_pi - DUt_pi_reverse_29ee8
   - IEX_pi[u]__TaxaB + IEX_pi[u]__TaxaB_reverse_aa8ed - IEX_pi[u]__TaxaC
   + IEX_pi[u]__TaxaC_reverse_acfe9 - IEX_pi[u]__TaxaD
   + IEX_pi[u]__TaxaD_reverse_925d1 = 0
 pi[fe]: - EX_pi[fe] + EX_pi[fe]_reverse_b0f41 + UFEt_pi
   - UFEt_pi_reverse_d3357 = 0
 pi[d]: - Diet_EX_pi[d] + Diet_EX_pi[d]_reverse_c8f98 - DUt_pi
   + DUt_pi_reverse_29ee8 = 0
 pyr[u]: - IEX_pyr[u]__TaxaA + IEX_pyr[u]__TaxaA_reverse_07e94 - UFEt_pyr
   + UFEt_pyr_reverse_0a6cc + DUt_pyr - DUt_pyr_reverse_30a17
   - IEX_pyr[u]__TaxaB + IEX_pyr[u]__TaxaB_reverse_bc30d
   - IEX_pyr[u]__TaxaC + IEX_pyr[u]__TaxaC_reverse_74aed
   - IEX_pyr[u]__TaxaD + IEX_pyr[u]__TaxaD_reverse_0e464 = 0
 pyr[fe]: - EX_pyr[fe] + EX_pyr[fe]_reverse_82c4f + UFEt_pyr
   - UFEt_pyr_reverse_0a6cc = 0
 pyr[d]: - Diet_EX_pyr[d] + Diet_EX_pyr[d]_reverse_9f7e9 - DUt_pyr
   + DUt_pyr_reverse_30a17 = 0
 succ[u]: - IEX_succ[u]__TaxaA + IEX_succ[u]__TaxaA_reverse_75565
   - UFEt_succ + UFEt_succ_reverse_86c10 + DUt_succ
   - DUt_succ_reverse_713b0 - IEX_succ[u]__TaxaB
   + IEX_succ[u]__TaxaB_reverse_119d9 - IEX_succ[u]__TaxaC
   + IEX_succ[u]__TaxaC_reverse_0ba26 - IEX_succ[u]__TaxaD
   + IEX_succ[u]__TaxaD_reverse_34a9d = 0
 succ[fe]: - EX_succ[fe] + EX_succ[fe]_reverse_356f5 + UFEt_succ
   - UFEt_succ_reverse_86c10 = 0
 succ[d]: - Diet_EX_succ[d] + Diet_EX_succ[d]_reverse_dd34a - DUt_succ
   + DUt_succ_reverse_713b0 = 0
 acald[c]__TaxaB: - ACALD__TaxaB + ACALD__TaxaB_reverse_0c342
   + ACALDt__TaxaB - ACALDt__TaxaB_reverse_8c029 + ALCD2x__TaxaB
   - ALCD2x__TaxaB_reverse_ade19 = 0
 coa[c]__TaxaB: - ACALD__TaxaB + ACALD__TaxaB_reverse_0c342 - AKGDH__TaxaB
   + AKGDH__TaxaB_reverse_3efa8 + 3.7478 biomass_ecoli__TaxaB
   - 3.7478 biomass_ecoli__TaxaB_reverse_0ec14 + CS__TaxaB
   - CS__TaxaB_reverse_f87a4 + MALS__TaxaB - MALS__TaxaB_reverse_4238a
   - PDH__TaxaB + PDH__TaxaB_reverse_a6e99 - PFL__TaxaB
   + PFL__TaxaB_reverse_04d04 + PTAr__TaxaB - PTAr__TaxaB_reverse_888fc
   - SUCOAS__TaxaB + SUCOAS__TaxaB_reverse_3f8e7 = 0
 nad[c]__TaxaB: - ACALD__TaxaB + ACALD__TaxaB_reverse_0c342 - AKGDH__TaxaB
   + AKGDH__TaxaB_reverse_3efa8 - ALCD2x__TaxaB
   + ALCD2x__TaxaB_reverse_ade19 - 3.547 biomass_ecoli__TaxaB
   + 3.547 biomass_ecoli__TaxaB_reverse_0ec14 - GAPD__TaxaB
   + GAPD__TaxaB_reverse_5d653 - LDH_D__TaxaB + LDH_D__TaxaB_reverse_a2c03
   - MDH__TaxaB + MDH__TaxaB_reverse_5ed0d - ME1__TaxaB
   + ME1__TaxaB_reverse_12347 + NADH16__TaxaB - NADH16__TaxaB_reverse_7b263
   - NADTRHD__TaxaB + NADTRHD__TaxaB_reverse_7e219 - PDH__TaxaB
   + PDH__TaxaB_reverse_a6e99 + THD2__TaxaB - THD2__TaxaB_reverse_c7723
   = 0
 accoa[c]__TaxaB: ACALD__TaxaB - ACALD__TaxaB_reverse_0c342
   - 3.7478 biomass_ecoli__TaxaB
   + 3.7478 biomass_ecoli__TaxaB_reverse_0ec14 - CS__TaxaB
   + CS__TaxaB_reverse_f87a4 - MALS__TaxaB + MALS__TaxaB_reverse_4238a
   + PDH__TaxaB - PDH__TaxaB_reverse_a6e99 + PFL__TaxaB
   - PFL__TaxaB_reverse_04d04 - PTAr__TaxaB + PTAr__TaxaB_reverse_888fc
   = 0
 h[c]__TaxaB: ACALD__TaxaB - ACALD__TaxaB_reverse_0c342 + ACt2r__TaxaB
   - ACt2r__TaxaB_reverse_7f31a + AKGt2r__TaxaB
   - AKGt2r__TaxaB_reverse_2ee07 + ALCD2x__TaxaB
   - ALCD2x__TaxaB_reverse_ade19 + ATPM__TaxaB - ATPM__TaxaB_reverse_fe91b
   + 3 ATPS4r__TaxaB - 3 ATPS4r__TaxaB_reverse_02fd2
   + 59.81 biomass_ecoli__TaxaB - 59.81 biomass_ecoli__TaxaB_reverse_0ec14
   + CS__TaxaB - CS__TaxaB_reverse_f87a4 - 2 CYTBD__TaxaB
   + 2 CYTBD__TaxaB_reverse_6f2d1 + D_LACt2__TaxaB
   - D_LACt2__TaxaB_reverse_e9195 + ETOHt2r__TaxaB
   - ETOHt2r__TaxaB_reverse_36ad9 + FORt2__TaxaB
   - FORt2__TaxaB_reverse_e4586 + 2 FUMt2_2__TaxaB
   - 2 FUMt2_2__TaxaB_reverse_3f590 + G6PDH2r__TaxaB
   - G6PDH2r__TaxaB_reverse_43665 + GAPD__TaxaB - GAPD__TaxaB_reverse_5d653
   + GLNS__TaxaB - GLNS__TaxaB_reverse_20b2d + GLNabc__TaxaB
   - GLNabc__TaxaB_reverse_ce96d + GLUDy__TaxaB
   - GLUDy__TaxaB_reverse_f03ea - GLUSy__TaxaB + GLUSy__TaxaB_reverse_ea33d
   + GLUt2r__TaxaB - GLUt2r__TaxaB_reverse_82bb4 + LDH_D__TaxaB
   - LDH_D__TaxaB_reverse_a2c03 + MALS__TaxaB - MALS__TaxaB_reverse_4238a
   + 2 MALt2_2__TaxaB - 2 MALt2_2__TaxaB_reverse_bdb3e + MDH__TaxaB
   - MDH__TaxaB_reverse_5ed0d - 4 NADH16__TaxaB
   + 4 NADH16__TaxaB_reverse_7b263 + PFK__TaxaB - PFK__TaxaB_reverse_9051a
   + PGL__TaxaB - PGL__TaxaB_reverse_bc58f + PIt2r__TaxaB
   - PIt2r__TaxaB_reverse_b31fe + PPC__TaxaB - PPC__TaxaB_reverse_2cfb9
   + 2 PPS__TaxaB - 2 PPS__TaxaB_reverse_d02ce - PYK__TaxaB
   + PYK__TaxaB_reverse_72465 + PYRt2__TaxaB - PYRt2__TaxaB_reverse_ed749
   + 2 SUCCt2_2__TaxaB - 2 SUCCt2_2__TaxaB_reverse_4bdaf + SUCCt3__TaxaB
   - SUCCt3__TaxaB_reverse_f38cd + 2 THD2__TaxaB
   - 2 THD2__TaxaB_reverse_c7723 = 0
 nadh[c]__TaxaB: ACALD__TaxaB - ACALD__TaxaB_reverse_0c342 + AKGDH__TaxaB
   - AKGDH__TaxaB_reverse_3efa8 + ALCD2x__TaxaB
   - ALCD2x__TaxaB_reverse_ade19 + 3.547 biomass_ecoli__TaxaB
   - 3.547 biomass_ecoli__TaxaB_reverse_0ec14 + GAPD__TaxaB
   - GAPD__TaxaB_reverse_5d653 + LDH_D__TaxaB - LDH_D__TaxaB_reverse_a2c03
   + MDH__TaxaB - MDH__TaxaB_reverse_5ed0d + ME1__TaxaB
   - ME1__TaxaB_reverse_12347 - NADH16__TaxaB + NADH16__TaxaB_reverse_7b263
   + NADTRHD__TaxaB - NADTRHD__TaxaB_reverse_7e219 + PDH__TaxaB
   - PDH__TaxaB_reverse_a6e99 - THD2__TaxaB + THD2__TaxaB_reverse_c7723
   = 0
 acald[u]__TaxaB: - ACALDt__TaxaB + ACALDt__TaxaB_reverse_8c029
   + IEX_acald[u]__TaxaB - IEX_acald[u]__TaxaB_reverse_2ec49 = 0
 ac[c]__TaxaB: - ACKr__TaxaB + ACKr__TaxaB_reverse_498be + ACt2r__TaxaB
   - ACt2r__TaxaB_reverse_7f31a = 0
 atp[c]__TaxaB: - ACKr__TaxaB + ACKr__TaxaB_reverse_498be - ADK1__TaxaB
   + ADK1__TaxaB_reverse_53760 - ATPM__TaxaB + ATPM__TaxaB_reverse_fe91b
   + ATPS4r__TaxaB - ATPS4r__TaxaB_reverse_02fd2
   - 59.81 biomass_ecoli__TaxaB + 59.81 biomass_ecoli__TaxaB_reverse_0ec14
   - GLNS__TaxaB + GLNS__TaxaB_reverse_20b2d - GLNabc__TaxaB
   + GLNabc__TaxaB_reverse_ce96d - PFK__TaxaB + PFK__TaxaB_reverse_9051a
   - PGK__TaxaB + PGK__TaxaB_reverse_e4ffb - PPCK__TaxaB
   + PPCK__TaxaB_reverse_086ef - PPS__TaxaB + PPS__TaxaB_reverse_d02ce
   + PYK__TaxaB - PYK__TaxaB_reverse_72465 - SUCOAS__TaxaB
   + SUCOAS__TaxaB_reverse_3f8e7 = 0
 actp[c]__TaxaB: ACKr__TaxaB - ACKr__TaxaB_reverse_498be + PTAr__TaxaB
   - PTAr__TaxaB_reverse_888fc = 0
 adp[c]__TaxaB: ACKr__TaxaB - ACKr__TaxaB_reverse_498be + 2 ADK1__TaxaB
   - 2 ADK1__TaxaB_reverse_53760 + ATPM__TaxaB - ATPM__TaxaB_reverse_fe91b
   - ATPS4r__TaxaB + ATPS4r__TaxaB_reverse_02fd2
   + 59.81 biomass_ecoli__TaxaB - 59.81 biomass_ecoli__TaxaB_reverse_0ec14
   + GLNS__TaxaB - GLNS__TaxaB_reverse_20b2d + GLNabc__TaxaB
   - GLNabc__TaxaB_reverse_ce96d + PFK__TaxaB - PFK__TaxaB_reverse_9051a
   + PGK__TaxaB - PGK__TaxaB_reverse_e4ffb + PPCK__TaxaB
   - PPCK__TaxaB_reverse_086ef - PYK__TaxaB + PYK__TaxaB_reverse_72465
   + SUCOAS__TaxaB - SUCOAS__TaxaB_reverse_3f8e7 = 0
 cit[c]__TaxaB: - ACONTa__TaxaB + ACONTa__TaxaB_reverse_7e72d + CS__TaxaB
   - CS__TaxaB_reverse_f87a4 = 0
 acon_C[c]__TaxaB: ACONTa__TaxaB - ACONTa__TaxaB_reverse_7e72d
   - ACONTb__TaxaB + ACONTb__TaxaB_reverse_93be1 = 0
 h2o[c]__TaxaB: ACONTa__TaxaB - ACONTa__TaxaB_reverse_7e72d - ACONTb__TaxaB
   + ACONTb__TaxaB_reverse_93be1 - ATPM__TaxaB + ATPM__TaxaB_reverse_fe91b
   + ATPS4r__TaxaB - ATPS4r__TaxaB_reverse_02fd2
   - 59.81 biomass_ecoli__TaxaB + 59.81 biomass_ecoli__TaxaB_reverse_0ec14
   - CS__TaxaB + CS__TaxaB_reverse_f87a4 + CYTBD__TaxaB
   - CYTBD__TaxaB_reverse_6f2d1 + ENO__TaxaB - ENO__TaxaB_reverse_e481f
   - FBP__TaxaB + FBP__TaxaB_reverse_bc783 - FUM__TaxaB
   + FUM__TaxaB_reverse_54296 - GLNabc__TaxaB + GLNabc__TaxaB_reverse_ce96d
   - GLUDy__TaxaB + GLUDy__TaxaB_reverse_f03ea - GLUN__TaxaB
   + GLUN__TaxaB_reverse_7f8f0 + H2Ot__TaxaB - H2Ot__TaxaB_reverse_b732d
   - MALS__TaxaB + MALS__TaxaB_reverse_4238a - PGL__TaxaB
   + PGL__TaxaB_reverse_bc58f - PPC__TaxaB + PPC__TaxaB_reverse_2cfb9
   - PPS__TaxaB + PPS__TaxaB_reverse_d02ce = 0
 icit[c]__TaxaB: ACONTb__TaxaB - ACONTb__TaxaB_reverse_93be1
   - ICDHyr__TaxaB + ICDHyr__TaxaB_reverse_fe281 - ICL__TaxaB
   + ICL__TaxaB_reverse_828bd = 0
 ac[u]__TaxaB: - ACt2r__TaxaB + ACt2r__TaxaB_reverse_7f31a
   + IEX_ac[u]__TaxaB - IEX_ac[u]__TaxaB_reverse_b168f = 0
 h[u]__TaxaB: - ACt2r__TaxaB + ACt2r__TaxaB_reverse_7f31a - AKGt2r__TaxaB
   + AKGt2r__TaxaB_reverse_2ee07 - 4 ATPS4r__TaxaB
   + 4 ATPS4r__TaxaB_reverse_02fd2 + 2 CYTBD__TaxaB
   - 2 CYTBD__TaxaB_reverse_6f2d1 - D_LACt2__TaxaB
   + D_LACt2__TaxaB_reverse_e9195 - ETOHt2r__TaxaB
   + ETOHt2r__TaxaB_reverse_36ad9 + IEX_h[u]__TaxaB
   - IEX_h[u]__TaxaB_reverse_a4431 - FORt2__TaxaB
   + FORt2__TaxaB_reverse_e4586 - 2 FUMt2_2__TaxaB
   + 2 FUMt2_2__TaxaB_reverse_3f590 - GLUt2r__TaxaB
   + GLUt2r__TaxaB_reverse_82bb4 - 2 MALt2_2__TaxaB
   + 2 MALt2_2__TaxaB_reverse_bdb3e + 3 NADH16__TaxaB
   - 3 NADH16__TaxaB_reverse_7b263 - PIt2r__TaxaB
   + PIt2r__TaxaB_reverse_b31fe - PYRt2__TaxaB + PYRt2__TaxaB_reverse_ed749
   - 2 SUCCt2_2__TaxaB + 2 SUCCt2_2__TaxaB_reverse_4bdaf - SUCCt3__TaxaB
   + SUCCt3__TaxaB_reverse_f38cd - 2 THD2__TaxaB
   + 2 THD2__TaxaB_reverse_c7723 = 0
 amp[c]__TaxaB: - ADK1__TaxaB + ADK1__TaxaB_reverse_53760 + PPS__TaxaB
   - PPS__TaxaB_reverse_d02ce = 0
 akg[c]__TaxaB: - AKGDH__TaxaB + AKGDH__TaxaB_reverse_3efa8 + AKGt2r__TaxaB
   - AKGt2r__TaxaB_reverse_2ee07 + 4.1182 biomass_ecoli__TaxaB
   - 4.1182 biomass_ecoli__TaxaB_reverse_0ec14 + GLUDy__TaxaB
   - GLUDy__TaxaB_reverse_f03ea - GLUSy__TaxaB + GLUSy__TaxaB_reverse_ea33d
   + ICDHyr__TaxaB - ICDHyr__TaxaB_reverse_fe281 = 0
 co2[c]__TaxaB: AKGDH__TaxaB - AKGDH__TaxaB_reverse_3efa8 + CO2t__TaxaB
   - CO2t__TaxaB_reverse_2f8f0 + GND__TaxaB - GND__TaxaB_reverse_d23bc
   + ICDHyr__TaxaB - ICDHyr__TaxaB_reverse_fe281 + ME1__TaxaB
   - ME1__TaxaB_reverse_12347 + ME2__TaxaB - ME2__TaxaB_reverse_63076
   + PDH__TaxaB - PDH__TaxaB_reverse_a6e99 - PPC__TaxaB
   + PPC__TaxaB_reverse_2cfb9 + PPCK__TaxaB - PPCK__TaxaB_reverse_086ef
   = 0
 succoa[c]__TaxaB: AKGDH__TaxaB - AKGDH__TaxaB_reverse_3efa8
   + SUCOAS__TaxaB - SUCOAS__TaxaB_reverse_3f8e7 = 0
 akg[u]__TaxaB: - AKGt2r__TaxaB + AKGt2r__TaxaB_reverse_2ee07
   + IEX_akg[u]__TaxaB - IEX_akg[u]__TaxaB_reverse_427b9 = 0
 etoh[c]__TaxaB: - ALCD2x__TaxaB + ALCD2x__TaxaB_reverse_ade19
   + ETOHt2r__TaxaB - ETOHt2r__TaxaB_reverse_36ad9 = 0
 pi[c]__TaxaB: ATPM__TaxaB - ATPM__TaxaB_reverse_fe91b - ATPS4r__TaxaB
   + ATPS4r__TaxaB_reverse_02fd2 + 59.81 biomass_ecoli__TaxaB
   - 59.81 biomass_ecoli__TaxaB_reverse_0ec14 + FBP__TaxaB
   - FBP__TaxaB_reverse_bc783 - GAPD__TaxaB + GAPD__TaxaB_reverse_5d653
   + GLNS__TaxaB - GLNS__TaxaB_reverse_20b2d + GLNabc__TaxaB
   - GLNabc__TaxaB_reverse_ce96d + PIt2r__TaxaB
   - PIt2r__TaxaB_reverse_b31fe + PPC__TaxaB - PPC__TaxaB_reverse_2cfb9
   + PPS__TaxaB - PPS__TaxaB_reverse_d02ce - PTAr__TaxaB
   + PTAr__TaxaB_reverse_888fc + SUCOAS__TaxaB
   - SUCOAS__TaxaB_reverse_3f8e7 = 0
 3pg[c]__TaxaB: - 1.496 biomass_ecoli__TaxaB
   + 1.496 biomass_ecoli__TaxaB_reverse_0ec14 - PGK__TaxaB
   + PGK__TaxaB_reverse_e4ffb + PGM__TaxaB - PGM__TaxaB_reverse_81ae9 = 0
 e4p[c]__TaxaB: - 0.361 biomass_ecoli__TaxaB
   + 0.361 biomass_ecoli__TaxaB_reverse_0ec14 + TALA__TaxaB
   - TALA__TaxaB_reverse_050f7 - TKT2__TaxaB + TKT2__TaxaB_reverse_39177
   = 0
 f6p[c]__TaxaB: - 0.0709 biomass_ecoli__TaxaB
   + 0.0709 biomass_ecoli__TaxaB_reverse_0ec14 + FBP__TaxaB
   - FBP__TaxaB_reverse_bc783 + FRUpts2__TaxaB
   - FRUpts2__TaxaB_reverse_d87d8 - PFK__TaxaB + PFK__TaxaB_reverse_9051a
   + PGI__TaxaB - PGI__TaxaB_reverse_ee07f + TALA__TaxaB
   - TALA__TaxaB_reverse_050f7 + TKT2__TaxaB - TKT2__TaxaB_reverse_39177
   = 0
 g3p[c]__TaxaB: - 0.129 biomass_ecoli__TaxaB
   + 0.129 biomass_ecoli__TaxaB_reverse_0ec14 + FBA__TaxaB
   - FBA__TaxaB_reverse_96afb - GAPD__TaxaB + GAPD__TaxaB_reverse_5d653
   - TALA__TaxaB + TALA__TaxaB_reverse_050f7 + TKT1__TaxaB
   - TKT1__TaxaB_reverse_d1c70 + TKT2__TaxaB - TKT2__TaxaB_reverse_39177
   + TPI__TaxaB - TPI__TaxaB_reverse_b9846 = 0
 g6p[c]__TaxaB: - 0.205 biomass_ecoli__TaxaB
   + 0.205 biomass_ecoli__TaxaB_reverse_0ec14 - G6PDH2r__TaxaB
   + G6PDH2r__TaxaB_reverse_43665 + GLCpts__TaxaB
   - GLCpts__TaxaB_reverse_ab812 - PGI__TaxaB + PGI__TaxaB_reverse_ee07f
   = 0
 gln__L[c]__TaxaB: - 0.2557 biomass_ecoli__TaxaB
   + 0.2557 biomass_ecoli__TaxaB_reverse_0ec14 + GLNS__TaxaB
   - GLNS__TaxaB_reverse_20b2d + GLNabc__TaxaB
   - GLNabc__TaxaB_reverse_ce96d - GLUN__TaxaB + GLUN__TaxaB_reverse_7f8f0
   - GLUSy__TaxaB + GLUSy__TaxaB_reverse_ea33d = 0
 glu__L[c]__TaxaB: - 4.9414 biomass_ecoli__TaxaB
   + 4.9414 biomass_ecoli__TaxaB_reverse_0ec14 - GLNS__TaxaB
   + GLNS__TaxaB_reverse_20b2d - GLUDy__TaxaB + GLUDy__TaxaB_reverse_f03ea
   + GLUN__TaxaB - GLUN__TaxaB_reverse_7f8f0 + 2 GLUSy__TaxaB
   - 2 GLUSy__TaxaB_reverse_ea33d + GLUt2r__TaxaB
   - GLUt2r__TaxaB_reverse_82bb4 = 0
 nadph[c]__TaxaB: - 13.0279 biomass_ecoli__TaxaB
   + 13.0279 biomass_ecoli__TaxaB_reverse_0ec14 + G6PDH2r__TaxaB
   - G6PDH2r__TaxaB_reverse_43665 + GLUDy__TaxaB
   - GLUDy__TaxaB_reverse_f03ea - GLUSy__TaxaB + GLUSy__TaxaB_reverse_ea33d
   + GND__TaxaB - GND__TaxaB_reverse_d23bc + ICDHyr__TaxaB
   - ICDHyr__TaxaB_reverse_fe281 + ME2__TaxaB - ME2__TaxaB_reverse_63076
   - NADTRHD__TaxaB + NADTRHD__TaxaB_reverse_7e219 + THD2__TaxaB
   - THD2__TaxaB_reverse_c7723 = 0
 oaa[c]__TaxaB: - 1.7867 biomass_ecoli__TaxaB
   + 1.7867 biomass_ecoli__TaxaB_reverse_0ec14 - CS__TaxaB
   + CS__TaxaB_reverse_f87a4 + MDH__TaxaB - MDH__TaxaB_reverse_5ed0d
   + PPC__TaxaB - PPC__TaxaB_reverse_2cfb9 - PPCK__TaxaB
   + PPCK__TaxaB_reverse_086ef = 0
 pep[c]__TaxaB: - 0.5191 biomass_ecoli__TaxaB
   + 0.5191 biomass_ecoli__TaxaB_reverse_0ec14 + ENO__TaxaB
   - ENO__TaxaB_reverse_e481f - FRUpts2__TaxaB
   + FRUpts2__TaxaB_reverse_d87d8 - GLCpts__TaxaB
   + GLCpts__TaxaB_reverse_ab812 - PPC__TaxaB + PPC__TaxaB_reverse_2cfb9
   + PPCK__TaxaB - PPCK__TaxaB_reverse_086ef + PPS__TaxaB
   - PPS__TaxaB_reverse_d02ce - PYK__TaxaB + PYK__TaxaB_reverse_72465 = 0
 pyr[c]__TaxaB: - 2.8328 biomass_ecoli__TaxaB
   + 2.8328 biomass_ecoli__TaxaB_reverse_0ec14 + FRUpts2__TaxaB
   - FRUpts2__TaxaB_reverse_d87d8 + GLCpts__TaxaB
   - GLCpts__TaxaB_reverse_ab812 + LDH_D__TaxaB
   - LDH_D__TaxaB_reverse_a2c03 + ME1__TaxaB - ME1__TaxaB_reverse_12347
   + ME2__TaxaB - ME2__TaxaB_reverse_63076 - PDH__TaxaB
   + PDH__TaxaB_reverse_a6e99 - PFL__TaxaB + PFL__TaxaB_reverse_04d04
   - PPS__TaxaB + PPS__TaxaB_reverse_d02ce + PYK__TaxaB
   - PYK__TaxaB_reverse_72465 + PYRt2__TaxaB - PYRt2__TaxaB_reverse_ed749
   = 0
 r5p[c]__TaxaB: - 0.8977 biomass_ecoli__TaxaB
   + 0.8977 biomass_ecoli__TaxaB_reverse_0ec14 - RPI__TaxaB
   + RPI__TaxaB_reverse_1475f - TKT1__TaxaB + TKT1__TaxaB_reverse_d1c70
   = 0
 nadp[c]__TaxaB: 13.0279 biomass_ecoli__TaxaB
   - 13.0279 biomass_ecoli__TaxaB_reverse_0ec14 - G6PDH2r__TaxaB
   + G6PDH2r__TaxaB_reverse_43665 - GLUDy__TaxaB
   + GLUDy__TaxaB_reverse_f03ea + GLUSy__TaxaB - GLUSy__TaxaB_reverse_ea33d
   - GND__TaxaB + GND__TaxaB_reverse_d23bc - ICDHyr__TaxaB
   + ICDHyr__TaxaB_reverse_fe281 - ME2__TaxaB + ME2__TaxaB_reverse_63076
   + NADTRHD__TaxaB - NADTRHD__TaxaB_reverse_7e219 - THD2__TaxaB
   + THD2__TaxaB_reverse_c7723 = 0
 biomass[c]__TaxaB: biomass_ecoli__TaxaB
   - biomass_ecoli__TaxaB_reverse_0ec14 - 0.2 communityBiomass
   + 0.2 communityBiomass_reverse_5ad9d = 0
 co2[u]__TaxaB: - CO2t__TaxaB + CO2t__TaxaB_reverse_2f8f0
   + IEX_co2[u]__TaxaB - IEX_co2[u]__TaxaB_reverse_60e91 = 0
 o2[c]__TaxaB: - 0.5 CYTBD__TaxaB + 0.5 CYTBD__TaxaB_reverse_6f2d1
   + O2t__TaxaB - O2t__TaxaB_reverse_ade13 = 0
 q8h2[c]__TaxaB: - CYTBD__TaxaB + CYTBD__TaxaB_reverse_6f2d1 - FRD7__TaxaB
   + FRD7__TaxaB_reverse_c16e9 + NADH16__TaxaB
   - NADH16__TaxaB_reverse_7b263 + SUCDi__TaxaB
   - SUCDi__TaxaB_reverse_c5ffa = 0
 q8[c]__TaxaB: CYTBD__TaxaB - CYTBD__TaxaB_reverse_6f2d1 + FRD7__TaxaB
   - FRD7__TaxaB_reverse_c16e9 - NADH16__TaxaB
   + NADH16__TaxaB_reverse_7b263 - SUCDi__TaxaB
   + SUCDi__TaxaB_reverse_c5ffa = 0
 lac__D[u]__TaxaB: - D_LACt2__TaxaB + D_LACt2__TaxaB_reverse_e9195
   + IEX_lac__D[u]__TaxaB - IEX_lac__D[u]__TaxaB_reverse_03348 = 0
 lac__D[c]__TaxaB: D_LACt2__TaxaB - D_LACt2__TaxaB_reverse_e9195
   - LDH_D__TaxaB + LDH_D__TaxaB_reverse_a2c03 = 0
 2pg[c]__TaxaB: - ENO__TaxaB + ENO__TaxaB_reverse_e481f - PGM__TaxaB
   + PGM__TaxaB_reverse_81ae9 = 0
 etoh[u]__TaxaB: - ETOHt2r__TaxaB + ETOHt2r__TaxaB_reverse_36ad9
   + EX[u]toh[u]__TaxaB - EX[u]toh[u]__TaxaB_reverse_0c398 = 0
 for[u]__TaxaB: IEX_for[u]__TaxaB - IEX_for[u]__TaxaB_reverse_10439
   - FORt2__TaxaB + FORt2__TaxaB_reverse_e4586 + FORti__TaxaB
   - FORti__TaxaB_reverse_1b19a = 0
 fru[u]__TaxaB: IEX_fru[u]__TaxaB - IEX_fru[u]__TaxaB_reverse_3b7e3
   - FRUpts2__TaxaB + FRUpts2__TaxaB_reverse_d87d8 = 0
 fum[u]__TaxaB: IEX_fum[u]__TaxaB - IEX_fum[u]__TaxaB_reverse_14c2d
   - FUMt2_2__TaxaB + FUMt2_2__TaxaB_reverse_3f590 = 0
 glc__D[u]__TaxaB: IEX_glc__D[u]__TaxaB
   - IEX_glc__D[u]__TaxaB_reverse_37fbc - GLCpts__TaxaB
   + GLCpts__TaxaB_reverse_ab812 = 0
 gln__L[u]__TaxaB: IEX_gln__L[u]__TaxaB
   - IEX_gln__L[u]__TaxaB_reverse_eb617 - GLNabc__TaxaB
   + GLNabc__TaxaB_reverse_ce96d = 0
 glu__L[u]__TaxaB: IEX_glu__L[u]__TaxaB
   - IEX_glu__L[u]__TaxaB_reverse_afc30 - GLUt2r__TaxaB
   + GLUt2r__TaxaB_reverse_82bb4 = 0
 h2o[u]__TaxaB: IEX_h2o[u]__TaxaB - IEX_h2o[u]__TaxaB_reverse_32bd0
   - H2Ot__TaxaB + H2Ot__TaxaB_reverse_b732d = 0
 mal__L[u]__TaxaB: IEX_mal__L[u]__TaxaB
   - IEX_mal__L[u]__TaxaB_reverse_e69bf - MALt2_2__TaxaB
   + MALt2_2__TaxaB_reverse_bdb3e = 0
 nh4[u]__TaxaB: IEX_nh4[u]__TaxaB - IEX_nh4[u]__TaxaB_reverse_e819b
   - NH4t__TaxaB + NH4t__TaxaB_reverse_a1555 = 0
 o2[u]__TaxaB: IEX_o2[u]__TaxaB - IEX_o2[u]__TaxaB_reverse_6c7da
   - O2t__TaxaB + O2t__TaxaB_reverse_ade13 = 0
 pi[u]__TaxaB: IEX_pi[u]__TaxaB - IEX_pi[u]__TaxaB_reverse_aa8ed
   - PIt2r__TaxaB + PIt2r__TaxaB_reverse_b31fe = 0
 pyr[u]__TaxaB: IEX_pyr[u]__TaxaB - IEX_pyr[u]__TaxaB_reverse_bc30d
   - PYRt2__TaxaB + PYRt2__TaxaB_reverse_ed749 = 0
 succ[u]__TaxaB: IEX_succ[u]__TaxaB - IEX_succ[u]__TaxaB_reverse_119d9
   - SUCCt2_2__TaxaB + SUCCt2_2__TaxaB_reverse_4bdaf + SUCCt3__TaxaB
   - SUCCt3__TaxaB_reverse_f38cd = 0
 fdp[c]__TaxaB: - FBA__TaxaB + FBA__TaxaB_reverse_96afb - FBP__TaxaB
   + FBP__TaxaB_reverse_bc783 + PFK__TaxaB - PFK__TaxaB_reverse_9051a = 0
 dhap[c]__TaxaB: FBA__TaxaB - FBA__TaxaB_reverse_96afb - TPI__TaxaB
   + TPI__TaxaB_reverse_b9846 = 0
 for[c]__TaxaB: FORt2__TaxaB - FORt2__TaxaB_reverse_e4586 - FORti__TaxaB
   + FORti__TaxaB_reverse_1b19a + PFL__TaxaB - PFL__TaxaB_reverse_04d04
   = 0
 fum[c]__TaxaB: - FRD7__TaxaB + FRD7__TaxaB_reverse_c16e9 - FUM__TaxaB
   + FUM__TaxaB_reverse_54296 + FUMt2_2__TaxaB
   - FUMt2_2__TaxaB_reverse_3f590 + SUCDi__TaxaB
   - SUCDi__TaxaB_reverse_c5ffa = 0
 succ[c]__TaxaB: FRD7__TaxaB - FRD7__TaxaB_reverse_c16e9 + ICL__TaxaB
   - ICL__TaxaB_reverse_828bd + SUCCt2_2__TaxaB
   - SUCCt2_2__TaxaB_reverse_4bdaf - SUCCt3__TaxaB
   + SUCCt3__TaxaB_reverse_f38cd - SUCDi__TaxaB
   + SUCDi__TaxaB_reverse_c5ffa - SUCOAS__TaxaB
   + SUCOAS__TaxaB_reverse_3f8e7 = 0
 mal__L[c]__TaxaB: FUM__TaxaB - FUM__TaxaB_reverse_54296 + MALS__TaxaB
   - MALS__TaxaB_reverse_4238a + MALt2_2__TaxaB
   - MALt2_2__TaxaB_reverse_bdb3e - MDH__TaxaB + MDH__TaxaB_reverse_5ed0d
   - ME1__TaxaB + ME1__TaxaB_reverse_12347 - ME2__TaxaB
   + ME2__TaxaB_reverse_63076 = 0
 6pgl[c]__TaxaB: G6PDH2r__TaxaB - G6PDH2r__TaxaB_reverse_43665 - PGL__TaxaB
   + PGL__TaxaB_reverse_bc58f = 0
 13dpg[c]__TaxaB: GAPD__TaxaB - GAPD__TaxaB_reverse_5d653 + PGK__TaxaB
   - PGK__TaxaB_reverse_e4ffb = 0
 nh4[c]__TaxaB: - GLNS__TaxaB + GLNS__TaxaB_reverse_20b2d + GLUDy__TaxaB
   - GLUDy__TaxaB_reverse_f03ea + GLUN__TaxaB - GLUN__TaxaB_reverse_7f8f0
   + NH4t__TaxaB - NH4t__TaxaB_reverse_a1555 = 0
 6pgc[c]__TaxaB: - GND__TaxaB + GND__TaxaB_reverse_d23bc + PGL__TaxaB
   - PGL__TaxaB_reverse_bc58f = 0
 ru5p__D[c]__TaxaB: GND__TaxaB - GND__TaxaB_reverse_d23bc - RPE__TaxaB
   + RPE__TaxaB_reverse_3c165 + RPI__TaxaB - RPI__TaxaB_reverse_1475f = 0
 glx[c]__TaxaB: ICL__TaxaB - ICL__TaxaB_reverse_828bd - MALS__TaxaB
   + MALS__TaxaB_reverse_4238a = 0
 xu5p__D[c]__TaxaB: RPE__TaxaB - RPE__TaxaB_reverse_3c165 - TKT1__TaxaB
   + TKT1__TaxaB_reverse_d1c70 - TKT2__TaxaB + TKT2__TaxaB_reverse_39177
   = 0
 s7p[c]__TaxaB: - TALA__TaxaB + TALA__TaxaB_reverse_050f7 + TKT1__TaxaB
   - TKT1__TaxaB_reverse_d1c70 = 0
 acald[c]__TaxaC: - ACALD__TaxaC + ACALD__TaxaC_reverse_b9484
   + ACALDt__TaxaC - ACALDt__TaxaC_reverse_d1aa6 + ALCD2x__TaxaC
   - ALCD2x__TaxaC_reverse_63b1c = 0
 coa[c]__TaxaC: - ACALD__TaxaC + ACALD__TaxaC_reverse_b9484 - AKGDH__TaxaC
   + AKGDH__TaxaC_reverse_a3c9d + 3.7478 biomass_ecoli__TaxaC
   - 3.7478 biomass_ecoli__TaxaC_reverse_41d18 + CS__TaxaC
   - CS__TaxaC_reverse_162ab + MALS__TaxaC - MALS__TaxaC_reverse_da103
   - PDH__TaxaC + PDH__TaxaC_reverse_87edc - PFL__TaxaC
   + PFL__TaxaC_reverse_3769b + PTAr__TaxaC - PTAr__TaxaC_reverse_48ac3
   - SUCOAS__TaxaC + SUCOAS__TaxaC_reverse_e07eb = 0
 nad[c]__TaxaC: - ACALD__TaxaC + ACALD__TaxaC_reverse_b9484 - AKGDH__TaxaC
   + AKGDH__TaxaC_reverse_a3c9d - ALCD2x__TaxaC
   + ALCD2x__TaxaC_reverse_63b1c - 3.547 biomass_ecoli__TaxaC
   + 3.547 biomass_ecoli__TaxaC_reverse_41d18 - GAPD__TaxaC
   + GAPD__TaxaC_reverse_17a1a - LDH_D__TaxaC + LDH_D__TaxaC_reverse_56e41
   - MDH__TaxaC + MDH__TaxaC_reverse_5bd24 - ME1__TaxaC
   + ME1__TaxaC_reverse_b9203 + NADH16__TaxaC - NADH16__TaxaC_reverse_fb2ba
   - NADTRHD__TaxaC + NADTRHD__TaxaC_reverse_ac758 - PDH__TaxaC
   + PDH__TaxaC_reverse_87edc + THD2__TaxaC - THD2__TaxaC_reverse_a7f74
   = 0
 accoa[c]__TaxaC: ACALD__TaxaC - ACALD__TaxaC_reverse_b9484
   - 3.7478 biomass_ecoli__TaxaC
   + 3.7478 biomass_ecoli__TaxaC_reverse_41d18 - CS__TaxaC
   + CS__TaxaC_reverse_162ab - MALS__TaxaC + MALS__TaxaC_reverse_da103
   + PDH__TaxaC - PDH__TaxaC_reverse_87edc + PFL__TaxaC
   - PFL__TaxaC_reverse_3769b - PTAr__TaxaC + PTAr__TaxaC_reverse_48ac3
   = 0
 h[c]__TaxaC: ACALD__TaxaC - ACALD__TaxaC_reverse_b9484 + ACt2r__TaxaC
   - ACt2r__TaxaC_reverse_b30b4 + AKGt2r__TaxaC
   - AKGt2r__TaxaC_reverse_72a4e + ALCD2x__TaxaC
   - ALCD2x__TaxaC_reverse_63b1c + ATPM__TaxaC - ATPM__TaxaC_reverse_d3438
   + 3 ATPS4r__TaxaC - 3 ATPS4r__TaxaC_reverse_edfcc
   + 59.81 biomass_ecoli__TaxaC - 59.81 biomass_ecoli__TaxaC_reverse_41d18
   + CS__TaxaC - CS__TaxaC_reverse_162ab - 2 CYTBD__TaxaC
   + 2 CYTBD__TaxaC_reverse_b5e26 + D_LACt2__TaxaC
   - D_LACt2__TaxaC_reverse_fb45c + ETOHt2r__TaxaC
   - ETOHt2r__TaxaC_reverse_5a87c + FORt2__TaxaC
   - FORt2__TaxaC_reverse_c2407 + 2 FUMt2_2__TaxaC
   - 2 FUMt2_2__TaxaC_reverse_4db1b + G6PDH2r__TaxaC
   - G6PDH2r__TaxaC_reverse_ab3e1 + GAPD__TaxaC - GAPD__TaxaC_reverse_17a1a
   + GLNS__TaxaC - GLNS__TaxaC_reverse_12565 + GLNabc__TaxaC
   - GLNabc__TaxaC_reverse_904d0 + GLUDy__TaxaC
   - GLUDy__TaxaC_reverse_990ca - GLUSy__TaxaC + GLUSy__TaxaC_reverse_84635
   + GLUt2r__TaxaC - GLUt2r__TaxaC_reverse_b3974 + LDH_D__TaxaC
   - LDH_D__TaxaC_reverse_56e41 + MALS__TaxaC - MALS__TaxaC_reverse_da103
   + 2 MALt2_2__TaxaC - 2 MALt2_2__TaxaC_reverse_680c0 + MDH__TaxaC
   - MDH__TaxaC_reverse_5bd24 - 4 NADH16__TaxaC
   + 4 NADH16__TaxaC_reverse_fb2ba + PFK__TaxaC - PFK__TaxaC_reverse_a987d
   + PGL__TaxaC - PGL__TaxaC_reverse_0181e + PIt2r__TaxaC
   - PIt2r__TaxaC_reverse_206b5 + PPC__TaxaC - PPC__TaxaC_reverse_b4752
   + 2 PPS__TaxaC - 2 PPS__TaxaC_reverse_dba7f - PYK__TaxaC
   + PYK__TaxaC_reverse_95865 + PYRt2__TaxaC - PYRt2__TaxaC_reverse_aadfb
   + 2 SUCCt2_2__TaxaC - 2 SUCCt2_2__TaxaC_reverse_572ac + SUCCt3__TaxaC
   - SUCCt3__TaxaC_reverse_a206f + 2 THD2__TaxaC
   - 2 THD2__TaxaC_reverse_a7f74 = 0
 nadh[c]__TaxaC: ACALD__TaxaC - ACALD__TaxaC_reverse_b9484 + AKGDH__TaxaC
   - AKGDH__TaxaC_reverse_a3c9d + ALCD2x__TaxaC
   - ALCD2x__TaxaC_reverse_63b1c + 3.547 biomass_ecoli__TaxaC
   - 3.547 biomass_ecoli__TaxaC_reverse_41d18 + GAPD__TaxaC
   - GAPD__TaxaC_reverse_17a1a + LDH_D__TaxaC - LDH_D__TaxaC_reverse_56e41
   + MDH__TaxaC - MDH__TaxaC_reverse_5bd24 + ME1__TaxaC
   - ME1__TaxaC_reverse_b9203 - NADH16__TaxaC + NADH16__TaxaC_reverse_fb2ba
   + NADTRHD__TaxaC - NADTRHD__TaxaC_reverse_ac758 + PDH__TaxaC
   - PDH__TaxaC_reverse_87edc - THD2__TaxaC + THD2__TaxaC_reverse_a7f74
   = 0
 acald[u]__TaxaC: - ACALDt__TaxaC + ACALDt__TaxaC_reverse_d1aa6
   + IEX_acald[u]__TaxaC - IEX_acald[u]__TaxaC_reverse_79433 = 0
 ac[c]__TaxaC: - ACKr__TaxaC + ACKr__TaxaC_reverse_be2cc + ACt2r__TaxaC
   - ACt2r__TaxaC_reverse_b30b4 = 0
 atp[c]__TaxaC: - ACKr__TaxaC + ACKr__TaxaC_reverse_be2cc - ADK1__TaxaC
   + ADK1__TaxaC_reverse_0a628 - ATPM__TaxaC + ATPM__TaxaC_reverse_d3438
   + ATPS4r__TaxaC - ATPS4r__TaxaC_reverse_edfcc
   - 59.81 biomass_ecoli__TaxaC + 59.81 biomass_ecoli__TaxaC_reverse_41d18
   - GLNS__TaxaC + GLNS__TaxaC_reverse_12565 - GLNabc__TaxaC
   + GLNabc__TaxaC_reverse_904d0 - PFK__TaxaC + PFK__TaxaC_reverse_a987d
   - PGK__TaxaC + PGK__TaxaC_reverse_aaf6a - PPCK__TaxaC
   + PPCK__TaxaC_reverse_400f8 - PPS__TaxaC + PPS__TaxaC_reverse_dba7f
   + PYK__TaxaC - PYK__TaxaC_reverse_95865 - SUCOAS__TaxaC
   + SUCOAS__TaxaC_reverse_e07eb = 0
 actp[c]__TaxaC: ACKr__TaxaC - ACKr__TaxaC_reverse_be2cc + PTAr__TaxaC
   - PTAr__TaxaC_reverse_48ac3 = 0
 adp[c]__TaxaC: ACKr__TaxaC - ACKr__TaxaC_reverse_be2cc + 2 ADK1__TaxaC
   - 2 ADK1__TaxaC_reverse_0a628 + ATPM__TaxaC - ATPM__TaxaC_reverse_d3438
   - ATPS4r__TaxaC + ATPS4r__TaxaC_reverse_edfcc
   + 59.81 biomass_ecoli__TaxaC - 59.81 biomass_ecoli__TaxaC_reverse_41d18
   + GLNS__TaxaC - GLNS__TaxaC_reverse_12565 + GLNabc__TaxaC
   - GLNabc__TaxaC_reverse_904d0 + PFK__TaxaC - PFK__TaxaC_reverse_a987d
   + PGK__TaxaC - PGK__TaxaC_reverse_aaf6a + PPCK__TaxaC
   - PPCK__TaxaC_reverse_400f8 - PYK__TaxaC + PYK__TaxaC_reverse_95865
   + SUCOAS__TaxaC - SUCOAS__TaxaC_reverse_e07eb = 0
 cit[c]__TaxaC: - ACONTa__TaxaC + ACONTa__TaxaC_reverse_21b27 + CS__TaxaC
   - CS__TaxaC_reverse_162ab = 0
 acon_C[c]__TaxaC: ACONTa__TaxaC - ACONTa__TaxaC_reverse_21b27
   - ACONTb__TaxaC + ACONTb__TaxaC_reverse_9f088 = 0
 h2o[c]__TaxaC: ACONTa__TaxaC - ACONTa__TaxaC_reverse_21b27 - ACONTb__TaxaC
   + ACONTb__TaxaC_reverse_9f088 - ATPM__TaxaC + ATPM__TaxaC_reverse_d3438
   + ATPS4r__TaxaC - ATPS4r__TaxaC_reverse_edfcc
   - 59.81 biomass_ecoli__TaxaC + 59.81 biomass_ecoli__TaxaC_reverse_41d18
   - CS__TaxaC + CS__TaxaC_reverse_162ab + CYTBD__TaxaC
   - CYTBD__TaxaC_reverse_b5e26 + ENO__TaxaC - ENO__TaxaC_reverse_e5c24
   - FBP__TaxaC + FBP__TaxaC_reverse_27fed - FUM__TaxaC
   + FUM__TaxaC_reverse_0d074 - GLNabc__TaxaC + GLNabc__TaxaC_reverse_904d0
   - GLUDy__TaxaC + GLUDy__TaxaC_reverse_990ca - GLUN__TaxaC
   + GLUN__TaxaC_reverse_d2cdc + H2Ot__TaxaC - H2Ot__TaxaC_reverse_3f623
   - MALS__TaxaC + MALS__TaxaC_reverse_da103 - PGL__TaxaC
   + PGL__TaxaC_reverse_0181e - PPC__TaxaC + PPC__TaxaC_reverse_b4752
   - PPS__TaxaC + PPS__TaxaC_reverse_dba7f = 0
 icit[c]__TaxaC: ACONTb__TaxaC - ACONTb__TaxaC_reverse_9f088
   - ICDHyr__TaxaC + ICDHyr__TaxaC_reverse_5612a - ICL__TaxaC
   + ICL__TaxaC_reverse_6b407 = 0
 ac[u]__TaxaC: - ACt2r__TaxaC + ACt2r__TaxaC_reverse_b30b4
   + IEX_ac[u]__TaxaC - IEX_ac[u]__TaxaC_reverse_82616 = 0
 h[u]__TaxaC: - ACt2r__TaxaC + ACt2r__TaxaC_reverse_b30b4 - AKGt2r__TaxaC
   + AKGt2r__TaxaC_reverse_72a4e - 4 ATPS4r__TaxaC
   + 4 ATPS4r__TaxaC_reverse_edfcc + 2 CYTBD__TaxaC
   - 2 CYTBD__TaxaC_reverse_b5e26 - D_LACt2__TaxaC
   + D_LACt2__TaxaC_reverse_fb45c - ETOHt2r__TaxaC
   + ETOHt2r__TaxaC_reverse_5a87c + IEX_h[u]__TaxaC
   - IEX_h[u]__TaxaC_reverse_b39fd - FORt2__TaxaC
   + FORt2__TaxaC_reverse_c2407 - 2 FUMt2_2__TaxaC
   + 2 FUMt2_2__TaxaC_reverse_4db1b - GLUt2r__TaxaC
   + GLUt2r__TaxaC_reverse_b3974 - 2 MALt2_2__TaxaC
   + 2 MALt2_2__TaxaC_reverse_680c0 + 3 NADH16__TaxaC
   - 3 NADH16__TaxaC_reverse_fb2ba - PIt2r__TaxaC
   + PIt2r__TaxaC_reverse_206b5 - PYRt2__TaxaC + PYRt2__TaxaC_reverse_aadfb
   - 2 SUCCt2_2__TaxaC + 2 SUCCt2_2__TaxaC_reverse_572ac - SUCCt3__TaxaC
   + SUCCt3__TaxaC_reverse_a206f - 2 THD2__TaxaC
   + 2 THD2__TaxaC_reverse_a7f74 = 0
 amp[c]__TaxaC: - ADK1__TaxaC + ADK1__TaxaC_reverse_0a628 + PPS__TaxaC
   - PPS__TaxaC_reverse_dba7f = 0
 akg[c]__TaxaC: - AKGDH__TaxaC + AKGDH__TaxaC_reverse_a3c9d + AKGt2r__TaxaC
   - AKGt2r__TaxaC_reverse_72a4e + 4.1182 biomass_ecoli__TaxaC
   - 4.1182 biomass_ecoli__TaxaC_reverse_41d18 + GLUDy__TaxaC
   - GLUDy__TaxaC_reverse_990ca - GLUSy__TaxaC + GLUSy__TaxaC_reverse_84635
   + ICDHyr__TaxaC - ICDHyr__TaxaC_reverse_5612a = 0
 co2[c]__TaxaC: AKGDH__TaxaC - AKGDH__TaxaC_reverse_a3c9d + CO2t__TaxaC
   - CO2t__TaxaC_reverse_3ea53 + GND__TaxaC - GND__TaxaC_reverse_87f3c
   + ICDHyr__TaxaC - ICDHyr__TaxaC_reverse_5612a + ME1__TaxaC
   - ME1__TaxaC_reverse_b9203 + ME2__TaxaC - ME2__TaxaC_reverse_cf032
   + PDH__TaxaC - PDH__TaxaC_reverse_87edc - PPC__TaxaC
   + PPC__TaxaC_reverse_b4752 + PPCK__TaxaC - PPCK__TaxaC_reverse_400f8
   = 0
 succoa[c]__TaxaC: AKGDH__TaxaC - AKGDH__TaxaC_reverse_a3c9d
   + SUCOAS__TaxaC - SUCOAS__TaxaC_reverse_e07eb = 0
 akg[u]__TaxaC: - AKGt2r__TaxaC + AKGt2r__TaxaC_reverse_72a4e
   + IEX_akg[u]__TaxaC - IEX_akg[u]__TaxaC_reverse_2d021 = 0
 etoh[c]__TaxaC: - ALCD2x__TaxaC + ALCD2x__TaxaC_reverse_63b1c
   + ETOHt2r__TaxaC - ETOHt2r__TaxaC_reverse_5a87c = 0
 pi[c]__TaxaC: ATPM__TaxaC - ATPM__TaxaC_reverse_d3438 - ATPS4r__TaxaC
   + ATPS4r__TaxaC_reverse_edfcc + 59.81 biomass_ecoli__TaxaC
   - 59.81 biomass_ecoli__TaxaC_reverse_41d18 + FBP__TaxaC
   - FBP__TaxaC_reverse_27fed - GAPD__TaxaC + GAPD__TaxaC_reverse_17a1a
   + GLNS__TaxaC - GLNS__TaxaC_reverse_12565 + GLNabc__TaxaC
   - GLNabc__TaxaC_reverse_904d0 + PIt2r__TaxaC
   - PIt2r__TaxaC_reverse_206b5 + PPC__TaxaC - PPC__TaxaC_reverse_b4752
   + PPS__TaxaC - PPS__TaxaC_reverse_dba7f - PTAr__TaxaC
   + PTAr__TaxaC_reverse_48ac3 + SUCOAS__TaxaC
   - SUCOAS__TaxaC_reverse_e07eb = 0
 3pg[c]__TaxaC: - 1.496 biomass_ecoli__TaxaC
   + 1.496 biomass_ecoli__TaxaC_reverse_41d18 - PGK__TaxaC
   + PGK__TaxaC_reverse_aaf6a + PGM__TaxaC - PGM__TaxaC_reverse_cfa90 = 0
 e4p[c]__TaxaC: - 0.361 biomass_ecoli__TaxaC
   + 0.361 biomass_ecoli__TaxaC_reverse_41d18 + TALA__TaxaC
   - TALA__TaxaC_reverse_8e9b1 - TKT2__TaxaC + TKT2__TaxaC_reverse_f5ba9
   = 0
 f6p[c]__TaxaC: - 0.0709 biomass_ecoli__TaxaC
   + 0.0709 biomass_ecoli__TaxaC_reverse_41d18 + FBP__TaxaC
   - FBP__TaxaC_reverse_27fed + FRUpts2__TaxaC
   - FRUpts2__TaxaC_reverse_01b00 - PFK__TaxaC + PFK__TaxaC_reverse_a987d
   + PGI__TaxaC - PGI__TaxaC_reverse_b32d8 + TALA__TaxaC
   - TALA__TaxaC_reverse_8e9b1 + TKT2__TaxaC - TKT2__TaxaC_reverse_f5ba9
   = 0
 g3p[c]__TaxaC: - 0.129 biomass_ecoli__TaxaC
   + 0.129 biomass_ecoli__TaxaC_reverse_41d18 + FBA__TaxaC
   - FBA__TaxaC_reverse_50114 - GAPD__TaxaC + GAPD__TaxaC_reverse_17a1a
   - TALA__TaxaC + TALA__TaxaC_reverse_8e9b1 + TKT1__TaxaC
   - TKT1__TaxaC_reverse_61d41 + TKT2__TaxaC - TKT2__TaxaC_reverse_f5ba9
   + TPI__TaxaC - TPI__TaxaC_reverse_6b386 = 0
 g6p[c]__TaxaC: - 0.205 biomass_ecoli__TaxaC
   + 0.205 biomass_ecoli__TaxaC_reverse_41d18 - G6PDH2r__TaxaC
   + G6PDH2r__TaxaC_reverse_ab3e1 + GLCpts__TaxaC
   - GLCpts__TaxaC_reverse_d4b50 - PGI__TaxaC + PGI__TaxaC_reverse_b32d8
   = 0
 gln__L[c]__TaxaC: - 0.2557 biomass_ecoli__TaxaC
   + 0.2557 biomass_ecoli__TaxaC_reverse_41d18 + GLNS__TaxaC
   - GLNS__TaxaC_reverse_12565 + GLNabc__TaxaC
   - GLNabc__TaxaC_reverse_904d0 - GLUN__TaxaC + GLUN__TaxaC_reverse_d2cdc
   - GLUSy__TaxaC + GLUSy__TaxaC_reverse_84635 = 0
 glu__L[c]__TaxaC: - 4.9414 biomass_ecoli__TaxaC
   + 4.9414 biomass_ecoli__TaxaC_reverse_41d18 - GLNS__TaxaC
   + GLNS__TaxaC_reverse_12565 - GLUDy__TaxaC + GLUDy__TaxaC_reverse_990ca
   + GLUN__TaxaC - GLUN__TaxaC_reverse_d2cdc + 2 GLUSy__TaxaC
   - 2 GLUSy__TaxaC_reverse_84635 + GLUt2r__TaxaC
   - GLUt2r__TaxaC_reverse_b3974 = 0
 nadph[c]__TaxaC: - 13.0279 biomass_ecoli__TaxaC
   + 13.0279 biomass_ecoli__TaxaC_reverse_41d18 + G6PDH2r__TaxaC
   - G6PDH2r__TaxaC_reverse_ab3e1 + GLUDy__TaxaC
   - GLUDy__TaxaC_reverse_990ca - GLUSy__TaxaC + GLUSy__TaxaC_reverse_84635
   + GND__TaxaC - GND__TaxaC_reverse_87f3c + ICDHyr__TaxaC
   - ICDHyr__TaxaC_reverse_5612a + ME2__TaxaC - ME2__TaxaC_reverse_cf032
   - NADTRHD__TaxaC + NADTRHD__TaxaC_reverse_ac758 + THD2__TaxaC
   - THD2__TaxaC_reverse_a7f74 = 0
 oaa[c]__TaxaC: - 1.7867 biomass_ecoli__TaxaC
   + 1.7867 biomass_ecoli__TaxaC_reverse_41d18 - CS__TaxaC
   + CS__TaxaC_reverse_162ab + MDH__TaxaC - MDH__TaxaC_reverse_5bd24
   + PPC__TaxaC - PPC__TaxaC_reverse_b4752 - PPCK__TaxaC
   + PPCK__TaxaC_reverse_400f8 = 0
 pep[c]__TaxaC: - 0.5191 biomass_ecoli__TaxaC
   + 0.5191 biomass_ecoli__TaxaC_reverse_41d18 + ENO__TaxaC
   - ENO__TaxaC_reverse_e5c24 - FRUpts2__TaxaC
   + FRUpts2__TaxaC_reverse_01b00 - GLCpts__TaxaC
   + GLCpts__TaxaC_reverse_d4b50 - PPC__TaxaC + PPC__TaxaC_reverse_b4752
   + PPCK__TaxaC - PPCK__TaxaC_reverse_400f8 + PPS__TaxaC
   - PPS__TaxaC_reverse_dba7f - PYK__TaxaC + PYK__TaxaC_reverse_95865 = 0
 pyr[c]__TaxaC: - 2.8328 biomass_ecoli__TaxaC
   + 2.8328 biomass_ecoli__TaxaC_reverse_41d18 + FRUpts2__TaxaC
   - FRUpts2__TaxaC_reverse_01b00 + GLCpts__TaxaC
   - GLCpts__TaxaC_reverse_d4b50 + LDH_D__TaxaC
   - LDH_D__TaxaC_reverse_56e41 + ME1__TaxaC - ME1__TaxaC_reverse_b9203
   + ME2__TaxaC - ME2__TaxaC_reverse_cf032 - PDH__TaxaC
   + PDH__TaxaC_reverse_87edc - PFL__TaxaC + PFL__TaxaC_reverse_3769b
   - PPS__TaxaC + PPS__TaxaC_reverse_dba7f + PYK__TaxaC
   - PYK__TaxaC_reverse_95865 + PYRt2__TaxaC - PYRt2__TaxaC_reverse_aadfb
   = 0
 r5p[c]__TaxaC: - 0.8977 biomass_ecoli__TaxaC
   + 0.8977 biomass_ecoli__TaxaC_reverse_41d18 - RPI__TaxaC
   + RPI__TaxaC_reverse_2ab80 - TKT1__TaxaC + TKT1__TaxaC_reverse_61d41
   = 0
 nadp[c]__TaxaC: 13.0279 biomass_ecoli__TaxaC
   - 13.0279 biomass_ecoli__TaxaC_reverse_41d18 - G6PDH2r__TaxaC
   + G6PDH2r__TaxaC_reverse_ab3e1 - GLUDy__TaxaC
   + GLUDy__TaxaC_reverse_990ca + GLUSy__TaxaC - GLUSy__TaxaC_reverse_84635
   - GND__TaxaC + GND__TaxaC_reverse_87f3c - ICDHyr__TaxaC
   + ICDHyr__TaxaC_reverse_5612a - ME2__TaxaC + ME2__TaxaC_reverse_cf032
   + NADTRHD__TaxaC - NADTRHD__TaxaC_reverse_ac758 - THD2__TaxaC
   + THD2__TaxaC_reverse_a7f74 = 0
 biomass[c]__TaxaC: biomass_ecoli__TaxaC
   - biomass_ecoli__TaxaC_reverse_41d18 - 0.3 communityBiomass
   + 0.3 communityBiomass_reverse_5ad9d = 0
 co2[u]__TaxaC: - CO2t__TaxaC + CO2t__TaxaC_reverse_3ea53
   + IEX_co2[u]__TaxaC - IEX_co2[u]__TaxaC_reverse_f8836 = 0
 o2[c]__TaxaC: - 0.5 CYTBD__TaxaC + 0.5 CYTBD__TaxaC_reverse_b5e26
   + O2t__TaxaC - O2t__TaxaC_reverse_6c905 = 0
 q8h2[c]__TaxaC: - CYTBD__TaxaC + CYTBD__TaxaC_reverse_b5e26 - FRD7__TaxaC
   + FRD7__TaxaC_reverse_eec45 + NADH16__TaxaC
   - NADH16__TaxaC_reverse_fb2ba + SUCDi__TaxaC
   - SUCDi__TaxaC_reverse_10c32 = 0
 q8[c]__TaxaC: CYTBD__TaxaC - CYTBD__TaxaC_reverse_b5e26 + FRD7__TaxaC
   - FRD7__TaxaC_reverse_eec45 - NADH16__TaxaC
   + NADH16__TaxaC_reverse_fb2ba - SUCDi__TaxaC
   + SUCDi__TaxaC_reverse_10c32 = 0
 lac__D[u]__TaxaC: - D_LACt2__TaxaC + D_LACt2__TaxaC_reverse_fb45c
   + IEX_lac__D[u]__TaxaC - IEX_lac__D[u]__TaxaC_reverse_baf61 = 0
 lac__D[c]__TaxaC: D_LACt2__TaxaC - D_LACt2__TaxaC_reverse_fb45c
   - LDH_D__TaxaC + LDH_D__TaxaC_reverse_56e41 = 0
 2pg[c]__TaxaC: - ENO__TaxaC + ENO__TaxaC_reverse_e5c24 - PGM__TaxaC
   + PGM__TaxaC_reverse_cfa90 = 0
 etoh[u]__TaxaC: - ETOHt2r__TaxaC + ETOHt2r__TaxaC_reverse_5a87c
   + EX[u]toh[u]__TaxaC - EX[u]toh[u]__TaxaC_reverse_65089 = 0
 for[u]__TaxaC: IEX_for[u]__TaxaC - IEX_for[u]__TaxaC_reverse_89537
   - FORt2__TaxaC + FORt2__TaxaC_reverse_c2407 + FORti__TaxaC
   - FORti__TaxaC_reverse_ef905 = 0
 fru[u]__TaxaC: IEX_fru[u]__TaxaC - IEX_fru[u]__TaxaC_reverse_d1617
   - FRUpts2__TaxaC + FRUpts2__TaxaC_reverse_01b00 = 0
 fum[u]__TaxaC: IEX_fum[u]__TaxaC - IEX_fum[u]__TaxaC_reverse_66a1c
   - FUMt2_2__TaxaC + FUMt2_2__TaxaC_reverse_4db1b = 0
 glc__D[u]__TaxaC: IEX_glc__D[u]__TaxaC
   - IEX_glc__D[u]__TaxaC_reverse_4dbc9 - GLCpts__TaxaC
   + GLCpts__TaxaC_reverse_d4b50 = 0
 gln__L[u]__TaxaC: IEX_gln__L[u]__TaxaC
   - IEX_gln__L[u]__TaxaC_reverse_d4f65 - GLNabc__TaxaC
   + GLNabc__TaxaC_reverse_904d0 = 0
 glu__L[u]__TaxaC: IEX_glu__L[u]__TaxaC
   - IEX_glu__L[u]__TaxaC_reverse_b053a - GLUt2r__TaxaC
   + GLUt2r__TaxaC_reverse_b3974 = 0
 h2o[u]__TaxaC: IEX_h2o[u]__TaxaC - IEX_h2o[u]__TaxaC_reverse_f3c2c
   - H2Ot__TaxaC + H2Ot__TaxaC_reverse_3f623 = 0
 mal__L[u]__TaxaC: IEX_mal__L[u]__TaxaC
   - IEX_mal__L[u]__TaxaC_reverse_2b6d1 - MALt2_2__TaxaC
   + MALt2_2__TaxaC_reverse_680c0 = 0
 nh4[u]__TaxaC: IEX_nh4[u]__TaxaC - IEX_nh4[u]__TaxaC_reverse_486f8
   - NH4t__TaxaC + NH4t__TaxaC_reverse_b46fe = 0
 o2[u]__TaxaC: IEX_o2[u]__TaxaC - IEX_o2[u]__TaxaC_reverse_9d59a
   - O2t__TaxaC + O2t__TaxaC_reverse_6c905 = 0
 pi[u]__TaxaC: IEX_pi[u]__TaxaC - IEX_pi[u]__TaxaC_reverse_acfe9
   - PIt2r__TaxaC + PIt2r__TaxaC_reverse_206b5 = 0
 pyr[u]__TaxaC: IEX_pyr[u]__TaxaC - IEX_pyr[u]__TaxaC_reverse_74aed
   - PYRt2__TaxaC + PYRt2__TaxaC_reverse_aadfb = 0
 succ[u]__TaxaC: IEX_succ[u]__TaxaC - IEX_succ[u]__TaxaC_reverse_0ba26
   - SUCCt2_2__TaxaC + SUCCt2_2__TaxaC_reverse_572ac + SUCCt3__TaxaC
   - SUCCt3__TaxaC_reverse_a206f = 0
 fdp[c]__TaxaC: - FBA__TaxaC + FBA__TaxaC_reverse_50114 - FBP__TaxaC
   + FBP__TaxaC_reverse_27fed + PFK__TaxaC - PFK__TaxaC_reverse_a987d = 0
 dhap[c]__TaxaC: FBA__TaxaC - FBA__TaxaC_reverse_50114 - TPI__TaxaC
   + TPI__TaxaC_reverse_6b386 = 0
 for[c]__TaxaC: FORt2__TaxaC - FORt2__TaxaC_reverse_c2407 - FORti__TaxaC
   + FORti__TaxaC_reverse_ef905 + PFL__TaxaC - PFL__TaxaC_reverse_3769b
   = 0
 fum[c]__TaxaC: - FRD7__TaxaC + FRD7__TaxaC_reverse_eec45 - FUM__TaxaC
   + FUM__TaxaC_reverse_0d074 + FUMt2_2__TaxaC
   - FUMt2_2__TaxaC_reverse_4db1b + SUCDi__TaxaC
   - SUCDi__TaxaC_reverse_10c32 = 0
 succ[c]__TaxaC: FRD7__TaxaC - FRD7__TaxaC_reverse_eec45 + ICL__TaxaC
   - ICL__TaxaC_reverse_6b407 + SUCCt2_2__TaxaC
   - SUCCt2_2__TaxaC_reverse_572ac - SUCCt3__TaxaC
   + SUCCt3__TaxaC_reverse_a206f - SUCDi__TaxaC
   + SUCDi__TaxaC_reverse_10c32 - SUCOAS__TaxaC
   + SUCOAS__TaxaC_reverse_e07eb = 0
 mal__L[c]__TaxaC: FUM__TaxaC - FUM__TaxaC_reverse_0d074 + MALS__TaxaC
   - MALS__TaxaC_reverse_da103 + MALt2_2__TaxaC
   - MALt2_2__TaxaC_reverse_680c0 - MDH__TaxaC + MDH__TaxaC_reverse_5bd24
   - ME1__TaxaC + ME1__TaxaC_reverse_b9203 - ME2__TaxaC
   + ME2__TaxaC_reverse_cf032 = 0
 6pgl[c]__TaxaC: G6PDH2r__TaxaC - G6PDH2r__TaxaC_reverse_ab3e1 - PGL__TaxaC
   + PGL__TaxaC_reverse_0181e = 0
 13dpg[c]__TaxaC: GAPD__TaxaC - GAPD__TaxaC_reverse_17a1a + PGK__TaxaC
   - PGK__TaxaC_reverse_aaf6a = 0
 nh4[c]__TaxaC: - GLNS__TaxaC + GLNS__TaxaC_reverse_12565 + GLUDy__TaxaC
   - GLUDy__TaxaC_reverse_990ca + GLUN__TaxaC - GLUN__TaxaC_reverse_d2cdc
   + NH4t__TaxaC - NH4t__TaxaC_reverse_b46fe = 0
 6pgc[c]__TaxaC: - GND__TaxaC + GND__TaxaC_reverse_87f3c + PGL__TaxaC
   - PGL__TaxaC_reverse_0181e = 0
 ru5p__D[c]__TaxaC: GND__TaxaC - GND__TaxaC_reverse_87f3c - RPE__TaxaC
   + RPE__TaxaC_reverse_8923a + RPI__TaxaC - RPI__TaxaC_reverse_2ab80 = 0
 glx[c]__TaxaC: ICL__TaxaC - ICL__TaxaC_reverse_6b407 - MALS__TaxaC
   + MALS__TaxaC_reverse_da103 = 0
 xu5p__D[c]__TaxaC: RPE__TaxaC - RPE__TaxaC_reverse_8923a - TKT1__TaxaC
   + TKT1__TaxaC_reverse_61d41 - TKT2__TaxaC + TKT2__TaxaC_reverse_f5ba9
   = 0
 s7p[c]__TaxaC: - TALA__TaxaC + TALA__TaxaC_reverse_8e9b1 + TKT1__TaxaC
   - TKT1__TaxaC_reverse_61d41 = 0
 acald[c]__TaxaD: - ACALD__TaxaD + ACALD__TaxaD_reverse_6d0be
   + ACALDt__TaxaD - ACALDt__TaxaD_reverse_0a331 + ALCD2x__TaxaD
   - ALCD2x__TaxaD_reverse_f3fcd = 0
 coa[c]__TaxaD: - ACALD__TaxaD + ACALD__TaxaD_reverse_6d0be - AKGDH__TaxaD
   + AKGDH__TaxaD_reverse_049a5 + 3.7478 biomass_ecoli__TaxaD
   - 3.7478 biomass_ecoli__TaxaD_reverse_71a87 + CS__TaxaD
   - CS__TaxaD_reverse_aaf14 + MALS__TaxaD - MALS__TaxaD_reverse_e8b65
   - PDH__TaxaD + PDH__TaxaD_reverse_a3fa6 - PFL__TaxaD
   + PFL__TaxaD_reverse_a5dcf + PTAr__TaxaD - PTAr__TaxaD_reverse_bb1d7
   - SUCOAS__TaxaD + SUCOAS__TaxaD_reverse_21a9f = 0
 nad[c]__TaxaD: - ACALD__TaxaD + ACALD__TaxaD_reverse_6d0be - AKGDH__TaxaD
   + AKGDH__TaxaD_reverse_049a5 - ALCD2x__TaxaD
   + ALCD2x__TaxaD_reverse_f3fcd - 3.547 biomass_ecoli__TaxaD
   + 3.547 biomass_ecoli__TaxaD_reverse_71a87 - GAPD__TaxaD
   + GAPD__TaxaD_reverse_20540 - LDH_D__TaxaD + LDH_D__TaxaD_reverse_49604
   - MDH__TaxaD + MDH__TaxaD_reverse_b8d94 - ME1__TaxaD
   + ME1__TaxaD_reverse_e55a1 + NADH16__TaxaD - NADH16__TaxaD_reverse_a5f47
   - NADTRHD__TaxaD + NADTRHD__TaxaD_reverse_22e5b - PDH__TaxaD
   + PDH__TaxaD_reverse_a3fa6 + THD2__TaxaD - THD2__TaxaD_reverse_7e9dc
   = 0
 accoa[c]__TaxaD: ACALD__TaxaD - ACALD__TaxaD_reverse_6d0be
   - 3.7478 biomass_ecoli__TaxaD
   + 3.7478 biomass_ecoli__TaxaD_reverse_71a87 - CS__TaxaD
   + CS__TaxaD_reverse_aaf14 - MALS__TaxaD + MALS__TaxaD_reverse_e8b65
   + PDH__TaxaD - PDH__TaxaD_reverse_a3fa6 + PFL__TaxaD
   - PFL__TaxaD_reverse_a5dcf - PTAr__TaxaD + PTAr__TaxaD_reverse_bb1d7
   = 0
 h[c]__TaxaD: ACALD__TaxaD - ACALD__TaxaD_reverse_6d0be + ACt2r__TaxaD
   - ACt2r__TaxaD_reverse_36c26 + AKGt2r__TaxaD
   - AKGt2r__TaxaD_reverse_2bddb + ALCD2x__TaxaD
   - ALCD2x__TaxaD_reverse_f3fcd + ATPM__TaxaD - ATPM__TaxaD_reverse_39785
   + 3 ATPS4r__TaxaD - 3 ATPS4r__TaxaD_reverse_5028a
   + 59.81 biomass_ecoli__TaxaD - 59.81 biomass_ecoli__TaxaD_reverse_71a87
   + CS__TaxaD - CS__TaxaD_reverse_aaf14 - 2 CYTBD__TaxaD
   + 2 CYTBD__TaxaD_reverse_4bb10 + D_LACt2__TaxaD
   - D_LACt2__TaxaD_reverse_b0332 + ETOHt2r__TaxaD
   - ETOHt2r__TaxaD_reverse_9310d + FORt2__TaxaD
   - FORt2__TaxaD_reverse_8e3d1 + 2 FUMt2_2__TaxaD
   - 2 FUMt2_2__TaxaD_reverse_d9367 + G6PDH2r__TaxaD
   - G6PDH2r__TaxaD_reverse_a3f2e + GAPD__TaxaD - GAPD__TaxaD_reverse_20540
   + GLNS__TaxaD - GLNS__TaxaD_reverse_5ed39 + GLNabc__TaxaD
   - GLNabc__TaxaD_reverse_85438 + GLUDy__TaxaD
   - GLUDy__TaxaD_reverse_203e2 - GLUSy__TaxaD + GLUSy__TaxaD_reverse_f2768
   + GLUt2r__TaxaD - GLUt2r__TaxaD_reverse_cedc6 + LDH_D__TaxaD
   - LDH_D__TaxaD_reverse_49604 + MALS__TaxaD - MALS__TaxaD_reverse_e8b65
   + 2 MALt2_2__TaxaD - 2 MALt2_2__TaxaD_reverse_cf2d6 + MDH__TaxaD
   - MDH__TaxaD_reverse_b8d94 - 4 NADH16__TaxaD
   + 4 NADH16__TaxaD_reverse_a5f47 + PFK__TaxaD - PFK__TaxaD_reverse_846b7
   + PGL__TaxaD - PGL__TaxaD_reverse_91fdc + PIt2r__TaxaD
   - PIt2r__TaxaD_reverse_93097 + PPC__TaxaD - PPC__TaxaD_reverse_b4029
   + 2 PPS__TaxaD - 2 PPS__TaxaD_reverse_fdd65 - PYK__TaxaD
   + PYK__TaxaD_reverse_70730 + PYRt2__TaxaD - PYRt2__TaxaD_reverse_20065
   + 2 SUCCt2_2__TaxaD - 2 SUCCt2_2__TaxaD_reverse_29f5f + SUCCt3__TaxaD
   - SUCCt3__TaxaD_reverse_282b4 + 2 THD2__TaxaD
   - 2 THD2__TaxaD_reverse_7e9dc = 0
 nadh[c]__TaxaD: ACALD__TaxaD - ACALD__TaxaD_reverse_6d0be + AKGDH__TaxaD
   - AKGDH__TaxaD_reverse_049a5 + ALCD2x__TaxaD
   - ALCD2x__TaxaD_reverse_f3fcd + 3.547 biomass_ecoli__TaxaD
   - 3.547 biomass_ecoli__TaxaD_reverse_71a87 + GAPD__TaxaD
   - GAPD__TaxaD_reverse_20540 + LDH_D__TaxaD - LDH_D__TaxaD_reverse_49604
   + MDH__TaxaD - MDH__TaxaD_reverse_b8d94 + ME1__TaxaD
   - ME1__TaxaD_reverse_e55a1 - NADH16__TaxaD + NADH16__TaxaD_reverse_a5f47
   + NADTRHD__TaxaD - NADTRHD__TaxaD_reverse_22e5b + PDH__TaxaD
   - PDH__TaxaD_reverse_a3fa6 - THD2__TaxaD + THD2__TaxaD_reverse_7e9dc
   = 0
 acald[u]__TaxaD: - ACALDt__TaxaD + ACALDt__TaxaD_reverse_0a331
   + IEX_acald[u]__TaxaD - IEX_acald[u]__TaxaD_reverse_1bc86 = 0
 ac[c]__TaxaD: - ACKr__TaxaD + ACKr__TaxaD_reverse_f63fc + ACt2r__TaxaD
   - ACt2r__TaxaD_reverse_36c26 = 0
 atp[c]__TaxaD: - ACKr__TaxaD + ACKr__TaxaD_reverse_f63fc - ADK1__TaxaD
   + ADK1__TaxaD_reverse_af0f9 - ATPM__TaxaD + ATPM__TaxaD_reverse_39785
   + ATPS4r__TaxaD - ATPS4r__TaxaD_reverse_5028a
   - 59.81 biomass_ecoli__TaxaD + 59.81 biomass_ecoli__TaxaD_reverse_71a87
   - GLNS__TaxaD + GLNS__TaxaD_reverse_5ed39 - GLNabc__TaxaD
   + GLNabc__TaxaD_reverse_85438 - PFK__TaxaD + PFK__TaxaD_reverse_846b7
   - PGK__TaxaD + PGK__TaxaD_reverse_15b8f - PPCK__TaxaD
   + PPCK__TaxaD_reverse_30ef8 - PPS__TaxaD + PPS__TaxaD_reverse_fdd65
   + PYK__TaxaD - PYK__TaxaD_reverse_70730 - SUCOAS__TaxaD
   + SUCOAS__TaxaD_reverse_21a9f = 0
 actp[c]__TaxaD: ACKr__TaxaD - ACKr__TaxaD_reverse_f63fc + PTAr__TaxaD
   - PTAr__TaxaD_reverse_bb1d7 = 0
 adp[c]__TaxaD: ACKr__TaxaD - ACKr__TaxaD_reverse_f63fc + 2 ADK1__TaxaD
   - 2 ADK1__TaxaD_reverse_af0f9 + ATPM__TaxaD - ATPM__TaxaD_reverse_39785
   - ATPS4r__TaxaD + ATPS4r__TaxaD_reverse_5028a
   + 59.81 biomass_ecoli__TaxaD - 59.81 biomass_ecoli__TaxaD_reverse_71a87
   + GLNS__TaxaD - GLNS__TaxaD_reverse_5ed39 + GLNabc__TaxaD
   - GLNabc__TaxaD_reverse_85438 + PFK__TaxaD - PFK__TaxaD_reverse_846b7
   + PGK__TaxaD - PGK__TaxaD_reverse_15b8f + PPCK__TaxaD
   - PPCK__TaxaD_reverse_30ef8 - PYK__TaxaD + PYK__TaxaD_reverse_70730
   + SUCOAS__TaxaD - SUCOAS__TaxaD_reverse_21a9f = 0
 cit[c]__TaxaD: - ACONTa__TaxaD + ACONTa__TaxaD_reverse_013c6 + CS__TaxaD
   - CS__TaxaD_reverse_aaf14 = 0
 acon_C[c]__TaxaD: ACONTa__TaxaD - ACONTa__TaxaD_reverse_013c6
   - ACONTb__TaxaD + ACONTb__TaxaD_reverse_799dc = 0
 h2o[c]__TaxaD: ACONTa__TaxaD - ACONTa__TaxaD_reverse_013c6 - ACONTb__TaxaD
   + ACONTb__TaxaD_reverse_799dc - ATPM__TaxaD + ATPM__TaxaD_reverse_39785
   + ATPS4r__TaxaD - ATPS4r__TaxaD_reverse_5028a
   - 59.81 biomass_ecoli__TaxaD + 59.81 biomass_ecoli__TaxaD_reverse_71a87
   - CS__TaxaD + CS__TaxaD_reverse_aaf14 + CYTBD__TaxaD
   - CYTBD__TaxaD_reverse_4bb10 + ENO__TaxaD - ENO__TaxaD_reverse_7d37f
   - FBP__TaxaD + FBP__TaxaD_reverse_5bb38 - FUM__TaxaD
   + FUM__TaxaD_reverse_3a7c0 - GLNabc__TaxaD + GLNabc__TaxaD_reverse_85438
   - GLUDy__TaxaD + GLUDy__TaxaD_reverse_203e2 - GLUN__TaxaD
   + GLUN__TaxaD_reverse_4524e + H2Ot__TaxaD - H2Ot__TaxaD_reverse_37165
   - MALS__TaxaD + MALS__TaxaD_reverse_e8b65 - PGL__TaxaD
   + PGL__TaxaD_reverse_91fdc - PPC__TaxaD + PPC__TaxaD_reverse_b4029
   - PPS__TaxaD + PPS__TaxaD_reverse_fdd65 = 0
 icit[c]__TaxaD: ACONTb__TaxaD - ACONTb__TaxaD_reverse_799dc
   - ICDHyr__TaxaD + ICDHyr__TaxaD_reverse_f264a - ICL__TaxaD
   + ICL__TaxaD_reverse_715d6 = 0
 ac[u]__TaxaD: - ACt2r__TaxaD + ACt2r__TaxaD_reverse_36c26
   + IEX_ac[u]__TaxaD - IEX_ac[u]__TaxaD_reverse_63ca5 = 0
 h[u]__TaxaD: - ACt2r__TaxaD + ACt2r__TaxaD_reverse_36c26 - AKGt2r__TaxaD
   + AKGt2r__TaxaD_reverse_2bddb - 4 ATPS4r__TaxaD
   + 4 ATPS4r__TaxaD_reverse_5028a + 2 CYTBD__TaxaD
   - 2 CYTBD__TaxaD_reverse_4bb10 - D_LACt2__TaxaD
   + D_LACt2__TaxaD_reverse_b0332 - ETOHt2r__TaxaD
   + ETOHt2r__TaxaD_reverse_9310d + IEX_h[u]__TaxaD
   - IEX_h[u]__TaxaD_reverse_3e6c7 - FORt2__TaxaD
   + FORt2__TaxaD_reverse_8e3d1 - 2 FUMt2_2__TaxaD
   + 2 FUMt2_2__TaxaD_reverse_d9367 - GLUt2r__TaxaD
   + GLUt2r__TaxaD_reverse_cedc6 - 2 MALt2_2__TaxaD
   + 2 MALt2_2__TaxaD_reverse_cf2d6 + 3 NADH16__TaxaD
   - 3 NADH16__TaxaD_reverse_a5f47 - PIt2r__TaxaD
   + PIt2r__TaxaD_reverse_93097 - PYRt2__TaxaD + PYRt2__TaxaD_reverse_20065
   - 2 SUCCt2_2__TaxaD + 2 SUCCt2_2__TaxaD_reverse_29f5f - SUCCt3__TaxaD
   + SUCCt3__TaxaD_reverse_282b4 - 2 THD2__TaxaD
   + 2 THD2__TaxaD_reverse_7e9dc = 0
 amp[c]__TaxaD: - ADK1__TaxaD + ADK1__TaxaD_reverse_af0f9 + PPS__TaxaD
   - PPS__TaxaD_reverse_fdd65 = 0
 akg[c]__TaxaD: - AKGDH__TaxaD + AKGDH__TaxaD_reverse_049a5 + AKGt2r__TaxaD
   - AKGt2r__TaxaD_reverse_2bddb + 4.1182 biomass_ecoli__TaxaD
   - 4.1182 biomass_ecoli__TaxaD_reverse_71a87 + GLUDy__TaxaD
   - GLUDy__TaxaD_reverse_203e2 - GLUSy__TaxaD + GLUSy__TaxaD_reverse_f2768
   + ICDHyr__TaxaD - ICDHyr__TaxaD_reverse_f264a = 0
 co2[c]__TaxaD: AKGDH__TaxaD - AKGDH__TaxaD_reverse_049a5 + CO2t__TaxaD
   - CO2t__TaxaD_reverse_f34ba + GND__TaxaD - GND__TaxaD_reverse_efd8f
   + ICDHyr__TaxaD - ICDHyr__TaxaD_reverse_f264a + ME1__TaxaD
   - ME1__TaxaD_reverse_e55a1 + ME2__TaxaD - ME2__TaxaD_reverse_d911c
   + PDH__TaxaD - PDH__TaxaD_reverse_a3fa6 - PPC__TaxaD
   + PPC__TaxaD_reverse_b4029 + PPCK__TaxaD - PPCK__TaxaD_reverse_30ef8
   = 0
 succoa[c]__TaxaD: AKGDH__TaxaD - AKGDH__TaxaD_reverse_049a5
   + SUCOAS__TaxaD - SUCOAS__TaxaD_reverse_21a9f = 0
 akg[u]__TaxaD: - AKGt2r__TaxaD + AKGt2r__TaxaD_reverse_2bddb
   + IEX_akg[u]__TaxaD - IEX_akg[u]__TaxaD_reverse_a37cb = 0
 etoh[c]__TaxaD: - ALCD2x__TaxaD + ALCD2x__TaxaD_reverse_f3fcd
   + ETOHt2r__TaxaD - ETOHt2r__TaxaD_reverse_9310d = 0
 pi[c]__TaxaD: ATPM__TaxaD - ATPM__TaxaD_reverse_39785 - ATPS4r__TaxaD
   + ATPS4r__TaxaD_reverse_5028a + 59.81 biomass_ecoli__TaxaD
   - 59.81 biomass_ecoli__TaxaD_reverse_71a87 + FBP__TaxaD
   - FBP__TaxaD_reverse_5bb38 - GAPD__TaxaD + GAPD__TaxaD_reverse_20540
   + GLNS__TaxaD - GLNS__TaxaD_reverse_5ed39 + GLNabc__TaxaD
   - GLNabc__TaxaD_reverse_85438 + PIt2r__TaxaD
   - PIt2r__TaxaD_reverse_93097 + PPC__TaxaD - PPC__TaxaD_reverse_b4029
   + PPS__TaxaD - PPS__TaxaD_reverse_fdd65 - PTAr__TaxaD
   + PTAr__TaxaD_reverse_bb1d7 + SUCOAS__TaxaD
   - SUCOAS__TaxaD_reverse_21a9f = 0
 3pg[c]__TaxaD: - 1.496 biomass_ecoli__TaxaD
   + 1.496 biomass_ecoli__TaxaD_reverse_71a87 - PGK__TaxaD
   + PGK__TaxaD_reverse_15b8f + PGM__TaxaD - PGM__TaxaD_reverse_732f2 = 0
 e4p[c]__TaxaD: - 0.361 biomass_ecoli__TaxaD
   + 0.361 biomass_ecoli__TaxaD_reverse_71a87 + TALA__TaxaD
   - TALA__TaxaD_reverse_d504b - TKT2__TaxaD + TKT2__TaxaD_reverse_68f4a
   = 0
 f6p[c]__TaxaD: - 0.0709 biomass_ecoli__TaxaD
   + 0.0709 biomass_ecoli__TaxaD_reverse_71a87 + FBP__TaxaD
   - FBP__TaxaD_reverse_5bb38 + FRUpts2__TaxaD
   - FRUpts2__TaxaD_reverse_61ab1 - PFK__TaxaD + PFK__TaxaD_reverse_846b7
   + PGI__TaxaD - PGI__TaxaD_reverse_6ff0e + TALA__TaxaD
   - TALA__TaxaD_reverse_d504b + TKT2__TaxaD - TKT2__TaxaD_reverse_68f4a
   = 0
 g3p[c]__TaxaD: - 0.129 biomass_ecoli__TaxaD
   + 0.129 biomass_ecoli__TaxaD_reverse_71a87 + FBA__TaxaD
   - FBA__TaxaD_reverse_4068a - GAPD__TaxaD + GAPD__TaxaD_reverse_20540
   - TALA__TaxaD + TALA__TaxaD_reverse_d504b + TKT1__TaxaD
   - TKT1__TaxaD_reverse_f64cc + TKT2__TaxaD - TKT2__TaxaD_reverse_68f4a
   + TPI__TaxaD - TPI__TaxaD_reverse_1efde = 0
 g6p[c]__TaxaD: - 0.205 biomass_ecoli__TaxaD
   + 0.205 biomass_ecoli__TaxaD_reverse_71a87 - G6PDH2r__TaxaD
   + G6PDH2r__TaxaD_reverse_a3f2e + GLCpts__TaxaD
   - GLCpts__TaxaD_reverse_e25bf - PGI__TaxaD + PGI__TaxaD_reverse_6ff0e
   = 0
 gln__L[c]__TaxaD: - 0.2557 biomass_ecoli__TaxaD
   + 0.2557 biomass_ecoli__TaxaD_reverse_71a87 + GLNS__TaxaD
   - GLNS__TaxaD_reverse_5ed39 + GLNabc__TaxaD
   - GLNabc__TaxaD_reverse_85438 - GLUN__TaxaD + GLUN__TaxaD_reverse_4524e
   - GLUSy__TaxaD + GLUSy__TaxaD_reverse_f2768 = 0
 glu__L[c]__TaxaD: - 4.9414 biomass_ecoli__TaxaD
   + 4.9414 biomass_ecoli__TaxaD_reverse_71a87 - GLNS__TaxaD
   + GLNS__TaxaD_reverse_5ed39 - GLUDy__TaxaD + GLUDy__TaxaD_reverse_203e2
   + GLUN__TaxaD - GLUN__TaxaD_reverse_4524e + 2 GLUSy__TaxaD
   - 2 GLUSy__TaxaD_reverse_f2768 + GLUt2r__TaxaD
   - GLUt2r__TaxaD_reverse_cedc6 = 0
 nadph[c]__TaxaD: - 13.0279 biomass_ecoli__TaxaD
   + 13.0279 biomass_ecoli__TaxaD_reverse_71a87 + G6PDH2r__TaxaD
   - G6PDH2r__TaxaD_reverse_a3f2e + GLUDy__TaxaD
   - GLUDy__TaxaD_reverse_203e2 - GLUSy__TaxaD + GLUSy__TaxaD_reverse_f2768
   + GND__TaxaD - GND__TaxaD_reverse_efd8f + ICDHyr__TaxaD
   - ICDHyr__TaxaD_reverse_f264a + ME2__TaxaD - ME2__TaxaD_reverse_d911c
   - NADTRHD__TaxaD + NADTRHD__TaxaD_reverse_22e5b + THD2__TaxaD
   - THD2__TaxaD_reverse_7e9dc = 0
 oaa[c]__TaxaD: - 1.7867 biomass_ecoli__TaxaD
   + 1.7867 biomass_ecoli__TaxaD_reverse_71a87 - CS__TaxaD
   + CS__TaxaD_reverse_aaf14 + MDH__TaxaD - MDH__TaxaD_reverse_b8d94
   + PPC__TaxaD - PPC__TaxaD_reverse_b4029 - PPCK__TaxaD
   + PPCK__TaxaD_reverse_30ef8 = 0
 pep[c]__TaxaD: - 0.5191 biomass_ecoli__TaxaD
   + 0.5191 biomass_ecoli__TaxaD_reverse_71a87 + ENO__TaxaD
   - ENO__TaxaD_reverse_7d37f - FRUpts2__TaxaD
   + FRUpts2__TaxaD_reverse_61ab1 - GLCpts__TaxaD
   + GLCpts__TaxaD_reverse_e25bf - PPC__TaxaD + PPC__TaxaD_reverse_b4029
   + PPCK__TaxaD - PPCK__TaxaD_reverse_30ef8 + PPS__TaxaD
   - PPS__TaxaD_reverse_fdd65 - PYK__TaxaD + PYK__TaxaD_reverse_70730 = 0
 pyr[c]__TaxaD: - 2.8328 biomass_ecoli__TaxaD
   + 2.8328 biomass_ecoli__TaxaD_reverse_71a87 + FRUpts2__TaxaD
   - FRUpts2__TaxaD_reverse_61ab1 + GLCpts__TaxaD
   - GLCpts__TaxaD_reverse_e25bf + LDH_D__TaxaD
   - LDH_D__TaxaD_reverse_49604 + ME1__TaxaD - ME1__TaxaD_reverse_e55a1
   + ME2__TaxaD - ME2__TaxaD_reverse_d911c - PDH__TaxaD
   + PDH__TaxaD_reverse_a3fa6 - PFL__TaxaD + PFL__TaxaD_reverse_a5dcf
   - PPS__TaxaD + PPS__TaxaD_reverse_fdd65 + PYK__TaxaD
   - PYK__TaxaD_reverse_70730 + PYRt2__TaxaD - PYRt2__TaxaD_reverse_20065
   = 0
 r5p[c]__TaxaD: - 0.8977 biomass_ecoli__TaxaD
   + 0.8977 biomass_ecoli__TaxaD_reverse_71a87 - RPI__TaxaD
   + RPI__TaxaD_reverse_d146e - TKT1__TaxaD + TKT1__TaxaD_reverse_f64cc
   = 0
 nadp[c]__TaxaD: 13.0279 biomass_ecoli__TaxaD
   - 13.0279 biomass_ecoli__TaxaD_reverse_71a87 - G6PDH2r__TaxaD
   + G6PDH2r__TaxaD_reverse_a3f2e - GLUDy__TaxaD
   + GLUDy__TaxaD_reverse_203e2 + GLUSy__TaxaD - GLUSy__TaxaD_reverse_f2768
   - GND__TaxaD + GND__TaxaD_reverse_efd8f - ICDHyr__TaxaD
   + ICDHyr__TaxaD_reverse_f264a - ME2__TaxaD + ME2__TaxaD_reverse_d911c
   + NADTRHD__TaxaD - NADTRHD__TaxaD_reverse_22e5b - THD2__TaxaD
   + THD2__TaxaD_reverse_7e9dc = 0
 biomass[c]__TaxaD: biomass_ecoli__TaxaD
   - biomass_ecoli__TaxaD_reverse_71a87 - 0.4 communityBiomass
   + 0.4 communityBiomass_reverse_5ad9d = 0
 co2[u]__TaxaD: - CO2t__TaxaD + CO2t__TaxaD_reverse_f34ba
   + IEX_co2[u]__TaxaD - IEX_co2[u]__TaxaD_reverse_ce75e = 0
 o2[c]__TaxaD: - 0.5 CYTBD__TaxaD + 0.5 CYTBD__TaxaD_reverse_4bb10
   + O2t__TaxaD - O2t__TaxaD_reverse_169fb = 0
 q8h2[c]__TaxaD: - CYTBD__TaxaD + CYTBD__TaxaD_reverse_4bb10 - FRD7__TaxaD
   + FRD7__TaxaD_reverse_2b38f + NADH16__TaxaD
   - NADH16__TaxaD_reverse_a5f47 + SUCDi__TaxaD
   - SUCDi__TaxaD_reverse_35a39 = 0
 q8[c]__TaxaD: CYTBD__TaxaD - CYTBD__TaxaD_reverse_4bb10 + FRD7__TaxaD
   - FRD7__TaxaD_reverse_2b38f - NADH16__TaxaD
   + NADH16__TaxaD_reverse_a5f47 - SUCDi__TaxaD
   + SUCDi__TaxaD_reverse_35a39 = 0
 lac__D[u]__TaxaD: - D_LACt2__TaxaD + D_LACt2__TaxaD_reverse_b0332
   + IEX_lac__D[u]__TaxaD - IEX_lac__D[u]__TaxaD_reverse_115c3 = 0
 lac__D[c]__TaxaD: D_LACt2__TaxaD - D_LACt2__TaxaD_reverse_b0332
   - LDH_D__TaxaD + LDH_D__TaxaD_reverse_49604 = 0
 2pg[c]__TaxaD: - ENO__TaxaD + ENO__TaxaD_reverse_7d37f - PGM__TaxaD
   + PGM__TaxaD_reverse_732f2 = 0
 etoh[u]__TaxaD: - ETOHt2r__TaxaD + ETOHt2r__TaxaD_reverse_9310d
   + EX[u]toh[u]__TaxaD - EX[u]toh[u]__TaxaD_reverse_1d69d = 0
 for[u]__TaxaD: IEX_for[u]__TaxaD - IEX_for[u]__TaxaD_reverse_0a866
   - FORt2__TaxaD + FORt2__TaxaD_reverse_8e3d1 + FORti__TaxaD
   - FORti__TaxaD_reverse_a472d = 0
 fru[u]__TaxaD: IEX_fru[u]__TaxaD - IEX_fru[u]__TaxaD_reverse_8c9bf
   - FRUpts2__TaxaD + FRUpts2__TaxaD_reverse_61ab1 = 0
 fum[u]__TaxaD: IEX_fum[u]__TaxaD - IEX_fum[u]__TaxaD_reverse_3415d
   - FUMt2_2__TaxaD + FUMt2_2__TaxaD_reverse_d9367 = 0
 glc__D[u]__TaxaD: IEX_glc__D[u]__TaxaD
   - IEX_glc__D[u]__TaxaD_reverse_405dc - GLCpts__TaxaD
   + GLCpts__TaxaD_reverse_e25bf = 0
 gln__L[u]__TaxaD: IEX_gln__L[u]__TaxaD
   - IEX_gln__L[u]__TaxaD_reverse_23702 - GLNabc__TaxaD
   + GLNabc__TaxaD_reverse_85438 = 0
 glu__L[u]__TaxaD: IEX_glu__L[u]__TaxaD
   - IEX_glu__L[u]__TaxaD_reverse_db289 - GLUt2r__TaxaD
   + GLUt2r__TaxaD_reverse_cedc6 = 0
 h2o[u]__TaxaD: IEX_h2o[u]__TaxaD - IEX_h2o[u]__TaxaD_reverse_c1c75
   - H2Ot__TaxaD + H2Ot__TaxaD_reverse_37165 = 0
 mal__L[u]__TaxaD: IEX_mal__L[u]__TaxaD
   - IEX_mal__L[u]__TaxaD_reverse_36f62 - MALt2_2__TaxaD
   + MALt2_2__TaxaD_reverse_cf2d6 = 0
 nh4[u]__TaxaD: IEX_nh4[u]__TaxaD - IEX_nh4[u]__TaxaD_reverse_a85a3
   - NH4t__TaxaD + NH4t__TaxaD_reverse_27ea9 = 0
 o2[u]__TaxaD: IEX_o2[u]__TaxaD - IEX_o2[u]__TaxaD_reverse_a1b64
   - O2t__TaxaD + O2t__TaxaD_reverse_169fb = 0
 pi[u]__TaxaD: IEX_pi[u]__TaxaD - IEX_pi[u]__TaxaD_reverse_925d1
   - PIt2r__TaxaD + PIt2r__TaxaD_reverse_93097 = 0
 pyr[u]__TaxaD: IEX_pyr[u]__TaxaD - IEX_pyr[u]__TaxaD_reverse_0e464
   - PYRt2__TaxaD + PYRt2__TaxaD_reverse_20065 = 0
 succ[u]__TaxaD: IEX_succ[u]__TaxaD - IEX_succ[u]__TaxaD_reverse_34a9d
   - SUCCt2_2__TaxaD + SUCCt2_2__TaxaD_reverse_29f5f + SUCCt3__TaxaD
   - SUCCt3__TaxaD_reverse_282b4 = 0
 fdp[c]__TaxaD: - FBA__TaxaD + FBA__TaxaD_reverse_4068a - FBP__TaxaD
   + FBP__TaxaD_reverse_5bb38 + PFK__TaxaD - PFK__TaxaD_reverse_846b7 = 0
 dhap[c]__TaxaD: FBA__TaxaD - FBA__TaxaD_reverse_4068a - TPI__TaxaD
   + TPI__TaxaD_reverse_1efde = 0
 for[c]__TaxaD: FORt2__TaxaD - FORt2__TaxaD_reverse_8e3d1 - FORti__TaxaD
   + FORti__TaxaD_reverse_a472d + PFL__TaxaD - PFL__TaxaD_reverse_a5dcf
   = 0
 fum[c]__TaxaD: - FRD7__TaxaD + FRD7__TaxaD_reverse_2b38f - FUM__TaxaD
   + FUM__TaxaD_reverse_3a7c0 + FUMt2_2__TaxaD
   - FUMt2_2__TaxaD_reverse_d9367 + SUCDi__TaxaD
   - SUCDi__TaxaD_reverse_35a39 = 0
 succ[c]__TaxaD: FRD7__TaxaD - FRD7__TaxaD_reverse_2b38f + ICL__TaxaD
   - ICL__TaxaD_reverse_715d6 + SUCCt2_2__TaxaD
   - SUCCt2_2__TaxaD_reverse_29f5f - SUCCt3__TaxaD
   + SUCCt3__TaxaD_reverse_282b4 - SUCDi__TaxaD
   + SUCDi__TaxaD_reverse_35a39 - SUCOAS__TaxaD
   + SUCOAS__TaxaD_reverse_21a9f = 0
 mal__L[c]__TaxaD: FUM__TaxaD - FUM__TaxaD_reverse_3a7c0 + MALS__TaxaD
   - MALS__TaxaD_reverse_e8b65 + MALt2_2__TaxaD
   - MALt2_2__TaxaD_reverse_cf2d6 - MDH__TaxaD + MDH__TaxaD_reverse_b8d94
   - ME1__TaxaD + ME1__TaxaD_reverse_e55a1 - ME2__TaxaD
   + ME2__TaxaD_reverse_d911c = 0
 6pgl[c]__TaxaD: G6PDH2r__TaxaD - G6PDH2r__TaxaD_reverse_a3f2e - PGL__TaxaD
   + PGL__TaxaD_reverse_91fdc = 0
 13dpg[c]__TaxaD: GAPD__TaxaD - GAPD__TaxaD_reverse_20540 + PGK__TaxaD
   - PGK__TaxaD_reverse_15b8f = 0
 nh4[c]__TaxaD: - GLNS__TaxaD + GLNS__TaxaD_reverse_5ed39 + GLUDy__TaxaD
   - GLUDy__TaxaD_reverse_203e2 + GLUN__TaxaD - GLUN__TaxaD_reverse_4524e
   + NH4t__TaxaD - NH4t__TaxaD_reverse_27ea9 = 0
 6pgc[c]__TaxaD: - GND__TaxaD + GND__TaxaD_reverse_efd8f + PGL__TaxaD
   - PGL__TaxaD_reverse_91fdc = 0
 ru5p__D[c]__TaxaD: GND__TaxaD - GND__TaxaD_reverse_efd8f - RPE__TaxaD
   + RPE__TaxaD_reverse_82d61 + RPI__TaxaD - RPI__TaxaD_reverse_d146e = 0
 glx[c]__TaxaD: ICL__TaxaD - ICL__TaxaD_reverse_715d6 - MALS__TaxaD
   + MALS__TaxaD_reverse_e8b65 = 0
 xu5p__D[c]__TaxaD: RPE__TaxaD - RPE__TaxaD_reverse_82d61 - TKT1__TaxaD
   + TKT1__TaxaD_reverse_f64cc - TKT2__TaxaD + TKT2__TaxaD_reverse_68f4a
   = 0
 s7p[c]__TaxaD: - TALA__TaxaD + TALA__TaxaD_reverse_d504b + TKT1__TaxaD
   - TKT1__TaxaD_reverse_f64cc = 0
 microbeBiomass[u]: - UFEt_microbeBiomass
   + UFEt_microbeBiomass_reverse_e9f5a + communityBiomass
   - communityBiomass_reverse_5ad9d = 0
 microbeBiomass[fe]: - EX_microbeBiomass[fe]
   + EX_microbeBiomass[fe]_reverse_784f8 + UFEt_microbeBiomass
   - UFEt_microbeBiomass_reverse_e9f5a = 0
Bounds
 ACALD__TaxaA <= 1000
 ACALD__TaxaA_reverse_3bd4d <= 1000
 ACALDt__TaxaA <= 1000
 ACALDt__TaxaA_reverse_7366b <= 1000
 ACKr__TaxaA <= 1000
 ACKr__TaxaA_reverse_35541 <= 1000
 ACONTa__TaxaA <= 1000
 ACONTa__TaxaA_reverse_dd81a <= 1000
 ACONTb__TaxaA <= 1000
 ACONTb__TaxaA_reverse_015ce <= 1000
 ACt2r__TaxaA <= 1000
 ACt2r__TaxaA_reverse_5e5e9 <= 1000
 ADK1__TaxaA <= 1000
 ADK1__TaxaA_reverse_af867 <= 1000
 AKGDH__TaxaA <= 1000
 AKGDH__TaxaA_reverse_48019 = 0
 AKGt2r__TaxaA <= 1000
 AKGt2r__TaxaA_reverse_56a5c <= 1000
 ALCD2x__TaxaA <= 1000
 ALCD2x__TaxaA_reverse_2678e <= 1000
 8.39 <= ATPM__TaxaA <= 1000
 ATPM__TaxaA_reverse_d6dff = 0
 ATPS4r__TaxaA <= 1000
 ATPS4r__TaxaA_reverse_4a5d8 <= 1000
 biomass_ecoli__TaxaA <= 1000
 biomass_ecoli__TaxaA_reverse_ed191 = 0
 CO2t__TaxaA <= 1000
 CO2t__TaxaA_reverse_5f188 <= 1000
 CS__TaxaA <= 1000
 CS__TaxaA_reverse_4b77c = 0
 CYTBD__TaxaA <= 1000
 CYTBD__TaxaA_reverse_983ef = 0
 D_LACt2__TaxaA <= 1000
 D_LACt2__TaxaA_reverse_a9e88 <= 1000
 ENO__TaxaA <= 1000
 ENO__TaxaA_reverse_b3da1 <= 1000
 ETOHt2r__TaxaA <= 1000
 ETOHt2r__TaxaA_reverse_8ca8d <= 1000
 IEX_ac[u]__TaxaA <= 1000
 IEX_ac[u]__TaxaA_reverse_d85aa <= 1000
 IEX_acald[u]__TaxaA <= 1000
 IEX_acald[u]__TaxaA_reverse_47467 <= 1000
 IEX_akg[u]__TaxaA <= 1000
 IEX_akg[u]__TaxaA_reverse_21614 <= 1000
 IEX_co2[u]__TaxaA <= 1000
 IEX_co2[u]__TaxaA_reverse_3ae77 <= 1000
 EX[u]toh[u]__TaxaA <= 1000
 EX[u]toh[u]__TaxaA_reverse_c58a3 <= 1000
 IEX_for[u]__TaxaA <= 1000
 IEX_for[u]__TaxaA_reverse_f0688 <= 1000
 IEX_fru[u]__TaxaA <= 1000
 IEX_fru[u]__TaxaA_reverse_0ac79 <= 1000
 IEX_fum[u]__TaxaA <= 1000
 IEX_fum[u]__TaxaA_reverse_b81bc <= 1000
 IEX_glc__D[u]__TaxaA <= 1000
 IEX_glc__D[u]__TaxaA_reverse_e84a1 <= 1000
 IEX_gln__L[u]__TaxaA <= 1000
 IEX_gln__L[u]__TaxaA_reverse_afdbf <= 1000
 IEX_glu__L[u]__TaxaA <= 1000
 IEX_glu__L[u]__TaxaA_reverse_2873c <= 1000
 IEX_h[u]__TaxaA <= 1000
 IEX_h[u]__TaxaA_reverse_29989 <= 1000
 IEX_h2o[u]__TaxaA <= 1000
 IEX_h2o[u]__TaxaA_reverse_dd098 <= 1000
 IEX_lac__D[u]__TaxaA <= 1000
 IEX_lac__D[u]__TaxaA_reverse_d234a <= 1000
 IEX_mal__L[u]__TaxaA <= 1000
 IEX_mal__L[u]__TaxaA_reverse_f589a <= 1000
 IEX_nh4[u]__TaxaA <= 1000
 IEX_nh4[u]__TaxaA_reverse_d1e33 <= 1000
 IEX_o2[u]__TaxaA <= 1000
 IEX_o2[u]__TaxaA_reverse_cbfe4 <= 1000
 IEX_pi[u]__TaxaA <= 1000
 IEX_pi[u]__TaxaA_reverse_621ea <= 1000
 IEX_pyr[u]__TaxaA <= 1000
 IEX_pyr[u]__TaxaA_reverse_07e94 <= 1000
 IEX_succ[u]__TaxaA <= 1000
 IEX_succ[u]__TaxaA_reverse_75565 <= 1000
 FBA__TaxaA <= 1000
 FBA__TaxaA_reverse_7aab1 <= 1000
 FBP__TaxaA <= 1000
 FBP__TaxaA_reverse_8adae = 0
 FORt2__TaxaA <= 1000
 FORt2__TaxaA_reverse_7aa37 = 0
 FORti__TaxaA <= 1000
 FORti__TaxaA_reverse_2b70c = 0
 FRD7__TaxaA <= 1000
 FRD7__TaxaA_reverse_80d05 = 0
 FRUpts2__TaxaA <= 1000
 FRUpts2__TaxaA_reverse_30526 = 0
 FUM__TaxaA <= 1000
 FUM__TaxaA_reverse_9925b <= 1000
 FUMt2_2__TaxaA <= 1000
 FUMt2_2__TaxaA_reverse_ae205 = 0
 G6PDH2r__TaxaA <= 1000
 G6PDH2r__TaxaA_reverse_798d9 <= 1000
 GAPD__TaxaA <= 1000
 GAPD__TaxaA_reverse_b06a4 <= 1000
 GLCpts__TaxaA <= 1000
 GLCpts__TaxaA_reverse_d6498 = 0
 GLNS__TaxaA <= 1000
 GLNS__TaxaA_reverse_b05d2 = 0
 GLNabc__TaxaA <= 1000
 GLNabc__TaxaA_reverse_317f3 = 0
 GLUDy__TaxaA <= 1000
 GLUDy__TaxaA_reverse_fb35f <= 1000
 GLUN__TaxaA <= 1000
 GLUN__TaxaA_reverse_85542 = 0
 GLUSy__TaxaA <= 1000
 GLUSy__TaxaA_reverse_ec3ea = 0
 GLUt2r__TaxaA <= 1000
 GLUt2r__TaxaA_reverse_286c6 <= 1000
 GND__TaxaA <= 1000
 GND__TaxaA_reverse_e35fe = 0
 H2Ot__TaxaA <= 1000
 H2Ot__TaxaA_reverse_127eb <= 1000
 ICDHyr__TaxaA <= 1000
 ICDHyr__TaxaA_reverse_54fda <= 1000
 ICL__TaxaA <= 1000
 ICL__TaxaA_reverse_b0fb1 = 0
 LDH_D__TaxaA <= 1000
 LDH_D__TaxaA_reverse_381f1 <= 1000
 MALS__TaxaA <= 1000
 MALS__TaxaA_reverse_c82e6 = 0
 MALt2_2__TaxaA <= 1000
 MALt2_2__TaxaA_reverse_87f6b = 0
 MDH__TaxaA <= 1000
 MDH__TaxaA_reverse_20572 <= 1000
 ME1__TaxaA <= 1000
 ME1__TaxaA_reverse_d8c3e = 0
 ME2__TaxaA <= 1000
 ME2__TaxaA_reverse_b7ffb = 0
 NADH16__TaxaA <= 1000
 NADH16__TaxaA_reverse_81d53 = 0
 NADTRHD__TaxaA <= 1000
 NADTRHD__TaxaA_reverse_2a766 = 0
 NH4t__TaxaA <= 1000
 NH4t__TaxaA_reverse_7a63f <= 1000
 O2t__TaxaA <= 1000
 O2t__TaxaA_reverse_b8a83 <= 1000
 PDH__TaxaA <= 1000
 PDH__TaxaA_reverse_cbf23 = 0
 PFK__TaxaA <= 1000
 PFK__TaxaA_reverse_ae474 = 0
 PFL__TaxaA <= 1000
 PFL__TaxaA_reverse_ff8bb = 0
 PGI__TaxaA <= 1000
 PGI__TaxaA_reverse_c02d7 <= 1000
 PGK__TaxaA <= 1000
 PGK__TaxaA_reverse_3f06d <= 1000
 PGL__TaxaA <= 1000
 PGL__TaxaA_reverse_23cba = 0
 PGM__TaxaA <= 1000
 PGM__TaxaA_reverse_b2abe <= 1000
 PIt2r__TaxaA <= 1000
 PIt2r__TaxaA_reverse_e3f73 <= 1000
 PPC__TaxaA <= 1000
 PPC__TaxaA_reverse_39721 = 0
 PPCK__TaxaA <= 1000
 PPCK__TaxaA_reverse_8d101 = 0
 PPS__TaxaA <= 1000
 PPS__TaxaA_reverse_422f8 = 0
 PTAr__TaxaA <= 1000
 PTAr__TaxaA_reverse_9e0a6 <= 1000
 PYK__TaxaA <= 1000
 PYK__TaxaA_reverse_dca3f = 0
 PYRt2__TaxaA <= 1000
 PYRt2__TaxaA_reverse_8f2e2 <= 1000
 RPE__TaxaA <= 1000
 RPE__TaxaA_reverse_50231 <= 1000
 RPI__TaxaA <= 1000
 RPI__TaxaA_reverse_b8683 <= 1000
 SUCCt2_2__TaxaA <= 1000
 SUCCt2_2__TaxaA_reverse_f20ba = 0
 SUCCt3__TaxaA <= 1000
 SUCCt3__TaxaA_reverse_f6546 = 0
 SUCDi__TaxaA <= 1000
 SUCDi__TaxaA_reverse_6fb91 = 0
 SUCOAS__TaxaA <= 1000
 SUCOAS__TaxaA_reverse_21627 <= 1000
 TALA__TaxaA <= 1000
 TALA__TaxaA_reverse_0bcac <= 1000
 THD2__TaxaA <= 1000
 THD2__TaxaA_reverse_f94fd = 0
 TKT1__TaxaA <= 1000
 TKT1__TaxaA_reverse_3bb31 <= 1000
 TKT2__TaxaA <= 1000
 TKT2__TaxaA_reverse_779d3 <= 1000
 TPI__TaxaA <= 1000
 TPI__TaxaA_reverse_9fc16 <= 1000
 EX_ac[fe] <= 1e+06
 EX_ac[fe]_reverse_15dd0 <= 1000
 UFEt_ac <= 1e+06
 UFEt_ac_reverse_0351f = 0
 Diet_EX_ac[d] <= 1000
 Diet_EX_ac[d]_reverse_c5b1c <= 1000
 DUt_ac <= 1e+06
 DUt_ac_reverse_d6d71 = 0
 EX_acald[fe] <= 1e+06
 EX_acald[fe]_reverse_745af <= 1000
 UFEt_acald <= 1e+06
 UFEt_acald_reverse_23948 = 0
 Diet_EX_acald[d] <= 1000
 Diet_EX_acald[d]_reverse_dcfe5 <= 1000
 DUt_acald <= 1e+06
 DUt_acald_reverse_88448 = 0
 EX_akg[fe] <= 1e+06
 EX_akg[fe]_reverse_3e084 <= 1000
 UFEt_akg <= 1e+06
 UFEt_akg_reverse_a80ce = 0
 Diet_EX_akg[d] <= 1000
 Diet_EX_akg[d]_reverse_e3523 <= 1000
 DUt_akg <= 1e+06
 DUt_akg_reverse_eb159 = 0
 EX_co2[fe] <= 1e+06
 EX_co2[fe]_reverse_9f050 <= 1000
 UFEt_co2 <= 1e+06
 UFEt_co2_reverse_de6bf = 0
 Diet_EX_co2[d] <= 1000
 Diet_EX_co2[d]_reverse_ec81a <= 1000
 DUt_co2 <= 1e+06
 DUt_co2_reverse_192a6 = 0
 EX_etoh[fe] <= 1e+06
 EX_etoh[fe]_reverse_406f6 <= 1000
 UFEt_etoh <= 1e+06
 UFEt_etoh_reverse_bbee7 = 0
 Diet_EX_etoh[d] <= 1000
 Diet_EX_etoh[d]_reverse_f4854 <= 1000
 DUt_etoh <= 1e+06
 DUt_etoh_reverse_c3f32 = 0
 EX_for[fe] <= 1e+06
 EX_for[fe]_reverse_5e89a <= 1000
 UFEt_for <= 1e+06
 UFEt_for_reverse_7662b = 0
 Diet_EX_for[d] <= 1000
 Diet_EX_for[d]_reverse_4b919 <= 1000
 DUt_for <= 1e+06
 DUt_for_reverse_f8c5e = 0
 EX_fru[fe] <= 1e+06
 EX_fru[fe]_reverse_9e214 <= 1000
 UFEt_fru <= 1e+06
 UFEt_fru_reverse_1d253 = 0
 Diet_EX_fru[d] <= 1000
 Diet_EX_fru[d]_reverse_079dc <= 1000
 DUt_fru <= 1e+06
 DUt_fru_reverse_b898a = 0
 EX_fum[fe] <= 1e+06
 EX_fum[fe]_reverse_7136c <= 1000
 UFEt_fum <= 1e+06
 UFEt_fum_reverse_9d01c = 0
 Diet_EX_fum[d] <= 1000
 Diet_EX_fum[d]_reverse_818ef <= 1000
 DUt_fum <= 1e+06
 DUt_fum_reverse_89aaf = 0
 EX_glc__D[fe] <= 1e+06
 EX_glc__D[fe]_reverse_48473 <= 1000
 UFEt_glc__D <= 1e+06
 UFEt_glc__D_reverse_691ca = 0
 Diet_EX_glc__D[d] <= 1000
 Diet_EX_glc__D[d]_reverse_a4b1e <= 1000
 DUt_glc__D <= 1e+06
 DUt_glc__D_reverse_c951b = 0
 EX_gln__L[fe] <= 1e+06
 EX_gln__L[fe]_reverse_f5e4f <= 1000
 UFEt_gln__L <= 1e+06
 UFEt_gln__L_reverse_76d9b = 0
 Diet_EX_gln__L[d] <= 1000
 Diet_EX_gln__L[d]_reverse_2ef88 <= 1000
 DUt_gln__L <= 1e+06
 DUt_gln__L_reverse_374f9 = 0
 EX_glu__L[fe] <= 1e+06
 EX_glu__L[fe]_reverse_24640 <= 1000
 UFEt_glu__L <= 1e+06
 UFEt_glu__L_reverse_e0c79 = 0
 Diet_EX_glu__L[d] <= 1000
 Diet_EX_glu__L[d]_reverse_6886b <= 1000
 DUt_glu__L <= 1e+06
 DUt_glu__L_reverse_331f7 = 0
 EX_h[fe] <= 1e+06
 EX_h[fe]_reverse_f5bf8 <= 1000
 UFEt_h <= 1e+06
 UFEt_h_reverse_0e86f = 0
 Diet_EX_h[d] <= 1000
 Diet_EX_h[d]_reverse_04d8d <= 1000
 DUt_h <= 1e+06
 DUt_h_reverse_e9a4a = 0
 EX_h2o[fe] <= 1e+06
 EX_h2o[fe]_reverse_e4723 <= 1000
 UFEt_h2o <= 1e+06
 UFEt_h2o_reverse_9fc26 = 0
 Diet_EX_h2o[d] <= 1000
 Diet_EX_h2o[d]_reverse_1812e <= 1000
 DUt_h2o <= 1e+06
 DUt_h2o_reverse_067ac = 0
 EX_lac__D[fe] <= 1e+06
 EX_lac__D[fe]_reverse_57587 <= 1000
 UFEt_lac__D <= 1e+06
 UFEt_lac__D_reverse_7be71 = 0
 Diet_EX_lac__D[d] <= 1000
 Diet_EX_lac__D[d]_reverse_07348 <= 1000
 DUt_lac__D <= 1e+06
 DUt_lac__D_reverse_cb570 = 0
 EX_mal__L[fe] <= 1e+06
 EX_mal__L[fe]_reverse_82659 <= 1000
 UFEt_mal__L <= 1e+06
 UFEt_mal__L_reverse_89431 = 0
 Diet_EX_mal__L[d] <= 1000
 Diet_EX_mal__L[d]_reverse_44e05 <= 1000
 DUt_mal__L <= 1e+06
 DUt_mal__L_reverse_4cd47 = 0
 EX_nh4[fe] <= 1e+06
 EX_nh4[fe]_reverse_02041 <= 1000
 UFEt_nh4 <= 1e+06
 UFEt_nh4_reverse_f9f09 = 0
 Diet_EX_nh4[d] <= 1000
 Diet_EX_nh4[d]_reverse_4a802 <= 1000
 DUt_nh4 <= 1e+06
 DUt_nh4_reverse_1b2f2 = 0
 EX_o2[fe] <= 1e+06
 EX_o2[fe]_reverse_4a084 <= 1000
 UFEt_o2 <= 1e+06
 UFEt_o2_reverse_d47cc = 0
 Diet_EX_o2[d] <= 1000
 Diet_EX_o2[d]_reverse_6f9b0 <= 1000
 DUt_o2 <= 1e+06
 DUt_o2_reverse_5b3e1 = 0
 EX_pi[fe] <= 1e+06
 EX_pi[fe]_reverse_b0f41 <= 1000
 UFEt_pi <= 1e+06
 UFEt_pi_reverse_d3357 = 0
 Diet_EX_pi[d] <= 1000
 Diet_EX_pi[d]_reverse_c8f98 <= 1000
 DUt_pi <= 1e+06
 DUt_pi_reverse_29ee8 = 0
 EX_pyr[fe] <= 1e+06
 EX_pyr[fe]_reverse_82c4f <= 1000
 UFEt_pyr <= 1e+06
 UFEt_pyr_reverse_0a6cc = 0
 Diet_EX_pyr[d] <= 1000
 Diet_EX_pyr[d]_reverse_9f7e9 <= 1000
 DUt_pyr <= 1e+06
 DUt_pyr_reverse_30a17 = 0
 EX_succ[fe] <= 1e+06
 EX_succ[fe]_reverse_356f5 <= 1000
 UFEt_succ <= 1e+06
 UFEt_succ_reverse_86c10 = 0
 Diet_EX_succ[d] <= 1000
 Diet_EX_succ[d]_reverse_dd34a <= 1000
 DUt_succ <= 1e+06
 DUt_succ_reverse_713b0 = 0
 ACALD__TaxaB <= 1000
 ACALD__TaxaB_reverse_0c342 <= 1000
 ACALDt__TaxaB <= 1000
 ACALDt__TaxaB_reverse_8c029 <= 1000
 ACKr__TaxaB <= 1000
 ACKr__TaxaB_reverse_498be <= 1000
 ACONTa__TaxaB <= 1000
 ACONTa__TaxaB_reverse_7e72d <= 1000
 ACONTb__TaxaB <= 1000
 ACONTb__TaxaB_reverse_93be1 <= 1000
 ACt2r__TaxaB <= 1000
 ACt2r__TaxaB_reverse_7f31a <= 1000
 ADK1__TaxaB <= 1000
 ADK1__TaxaB_reverse_53760 <= 1000
 AKGDH__TaxaB <= 1000
 AKGDH__TaxaB_reverse_3efa8 = 0
 AKGt2r__TaxaB <= 1000
 AKGt2r__TaxaB_reverse_2ee07 <= 1000
 ALCD2x__TaxaB <= 1000
 ALCD2x__TaxaB_reverse_ade19 <= 1000
 8.39 <= ATPM__TaxaB <= 1000
 ATPM__TaxaB_reverse_fe91b = 0
 ATPS4r__TaxaB <= 1000
 ATPS4r__TaxaB_reverse_02fd2 <= 1000
 biomass_ecoli__TaxaB <= 1000
 biomass_ecoli__TaxaB_reverse_0ec14 = 0
 CO2t__TaxaB <= 1000
 CO2t__TaxaB_reverse_2f8f0 <= 1000
 CS__TaxaB <= 1000
 CS__TaxaB_reverse_f87a4 = 0
 CYTBD__TaxaB <= 1000
 CYTBD__TaxaB_reverse_6f2d1 = 0
 D_LACt2__TaxaB <= 1000
 D_LACt2__TaxaB_reverse_e9195 <= 1000
 ENO__TaxaB <= 1000
 ENO__TaxaB_reverse_e481f <= 1000
 ETOHt2r__TaxaB <= 1000
 ETOHt2r__TaxaB_reverse_36ad9 <= 1000
 IEX_ac[u]__TaxaB <= 1000
 IEX_ac[u]__TaxaB_reverse_b168f <= 1000
 IEX_acald[u]__TaxaB <= 1000
 IEX_acald[u]__TaxaB_reverse_2ec49 <= 1000
 IEX_akg[u]__TaxaB <= 1000
 IEX_akg[u]__TaxaB_reverse_427b9 <= 1000
 IEX_co2[u]__TaxaB <= 1000
 IEX_co2[u]__TaxaB_reverse_60e91 <= 1000
 EX[u]toh[u]__TaxaB <= 1000
 EX[u]toh[u]__TaxaB_reverse_0c398 <= 1000
 IEX_for[u]__TaxaB <= 1000
 IEX_for[u]__TaxaB_reverse_10439 <= 1000
 IEX_fru[u]__TaxaB <= 1000
 IEX_fru[u]__TaxaB_reverse_3b7e3 <= 1000
 IEX_fum[u]__TaxaB <= 1000
 IEX_fum[u]__TaxaB_reverse_14c2d <= 1000
 IEX_glc__D[u]__TaxaB <= 1000
 IEX_glc__D[u]__TaxaB_reverse_37fbc <= 1000
 IEX_gln__L[u]__TaxaB <= 1000
 IEX_gln__L[u]__TaxaB_reverse_eb617 <= 1000
 IEX_glu__L[u]__TaxaB <= 1000
 IEX_glu__L[u]__TaxaB_reverse_afc30 <= 1000
 IEX_h[u]__TaxaB <= 1000
 IEX_h[u]__TaxaB_reverse_a4431 <= 1000
 IEX_h2o[u]__TaxaB <= 1000
 IEX_h2o[u]__TaxaB_reverse_32bd0 <= 1000
 IEX_lac__D[u]__TaxaB <= 1000
 IEX_lac__D[u]__TaxaB_reverse_03348 <= 1000
 IEX_mal__L[u]__TaxaB <= 1000
 IEX_mal__L[u]__TaxaB_reverse_e69bf <= 1000
 IEX_nh4[u]__TaxaB <= 1000
 IEX_nh4[u]__TaxaB_reverse_e819b <= 1000
 IEX_o2[u]__TaxaB <= 1000
 IEX_o2[u]__TaxaB_reverse_6c7da <= 1000
 IEX_pi[u]__TaxaB <= 1000
 IEX_pi[u]__TaxaB_reverse_aa8ed <= 1000
 IEX_pyr[u]__TaxaB <= 1000
 IEX_pyr[u]__TaxaB_reverse_bc30d <= 1000
 IEX_succ[u]__TaxaB <= 1000
 IEX_succ[u]__TaxaB_reverse_119d9 <= 1000
 FBA__TaxaB <= 1000
 FBA__TaxaB_reverse_96afb <= 1000
 FBP__TaxaB <= 1000
 FBP__TaxaB_reverse_bc783 = 0
 FORt2__TaxaB <= 1000
 FORt2__TaxaB_reverse_e4586 = 0
 FORti__TaxaB <= 1000
 FORti__TaxaB_reverse_1b19a = 0
 FRD7__TaxaB <= 1000
 FRD7__TaxaB_reverse_c16e9 = 0
 FRUpts2__TaxaB <= 1000
 FRUpts2__TaxaB_reverse_d87d8 = 0
 FUM__TaxaB <= 1000
 FUM__TaxaB_reverse_54296 <= 1000
 FUMt2_2__TaxaB <= 1000
 FUMt2_2__TaxaB_reverse_3f590 = 0
 G6PDH2r__TaxaB <= 1000
 G6PDH2r__TaxaB_reverse_43665 <= 1000
 GAPD__TaxaB <= 1000
 GAPD__TaxaB_reverse_5d653 <= 1000
 GLCpts__TaxaB <= 1000
 GLCpts__TaxaB_reverse_ab812 = 0
 GLNS__TaxaB <= 1000
 GLNS__TaxaB_reverse_20b2d = 0
 GLNabc__TaxaB <= 1000
 GLNabc__TaxaB_reverse_ce96d = 0
 GLUDy__TaxaB <= 1000
 GLUDy__TaxaB_reverse_f03ea <= 1000
 GLUN__TaxaB <= 1000
 GLUN__TaxaB_reverse_7f8f0 = 0
 GLUSy__TaxaB <= 1000
 GLUSy__TaxaB_reverse_ea33d = 0
 GLUt2r__TaxaB <= 1000
 GLUt2r__TaxaB_reverse_82bb4 <= 1000
 GND__TaxaB <= 1000
 GND__TaxaB_reverse_d23bc = 0
 H2Ot__TaxaB <= 1000
 H2Ot__TaxaB_reverse_b732d <= 1000
 ICDHyr__TaxaB <= 1000
 ICDHyr__TaxaB_reverse_fe281 <= 1000
 ICL__TaxaB <= 1000
 ICL__TaxaB_reverse_828bd = 0
 LDH_D__TaxaB <= 1000
 LDH_D__TaxaB_reverse_a2c03 <= 1000
 MALS__TaxaB <= 1000
 MALS__TaxaB_reverse_4238a = 0
 MALt2_2__TaxaB <= 1000
 MALt2_2__TaxaB_reverse_bdb3e = 0
 MDH__TaxaB <= 1000
 MDH__TaxaB_reverse_5ed0d <= 1000
 ME1__TaxaB <= 1000
 ME1__TaxaB_reverse_12347 = 0
 ME2__TaxaB <= 1000
 ME2__TaxaB_reverse_63076 = 0
 NADH16__TaxaB <= 1000
 NADH16__TaxaB_reverse_7b263 = 0
 NADTRHD__TaxaB <= 1000
 NADTRHD__TaxaB_reverse_7e219 = 0
 NH4t__TaxaB <= 1000
 NH4t__TaxaB_reverse_a1555 <= 1000
 O2t__TaxaB <= 1000
 O2t__TaxaB_reverse_ade13 <= 1000
 PDH__TaxaB <= 1000
 PDH__TaxaB_reverse_a6e99 = 0
 PFK__TaxaB <= 1000
 PFK__TaxaB_reverse_9051a = 0
 PFL__TaxaB <= 1000
 PFL__TaxaB_reverse_04d04 = 0
 PGI__TaxaB <= 1000
 PGI__TaxaB_reverse_ee07f <= 1000
 PGK__TaxaB <= 1000
 PGK__TaxaB_reverse_e4ffb <= 1000
 PGL__TaxaB <= 1000
 PGL__TaxaB_reverse_bc58f = 0
 PGM__TaxaB <= 1000
 PGM__TaxaB_reverse_81ae9 <= 1000
 PIt2r__TaxaB <= 1000
 PIt2r__TaxaB_reverse_b31fe <= 1000
 PPC__TaxaB <= 1000
 PPC__TaxaB_reverse_2cfb9 = 0
 PPCK__TaxaB <= 1000
 PPCK__TaxaB_reverse_086ef = 0
 PPS__TaxaB <= 1000
 PPS__TaxaB_reverse_d02ce = 0
 PTAr__TaxaB <= 1000
 PTAr__TaxaB_reverse_888fc <= 1000
 PYK__TaxaB <= 1000
 PYK__TaxaB_reverse_72465 = 0
 PYRt2__TaxaB <= 1000
 PYRt2__TaxaB_reverse_ed749 <= 1000
 RPE__TaxaB <= 1000
 RPE__TaxaB_reverse_3c165 <= 1000
 RPI__TaxaB <= 1000
 RPI__TaxaB_reverse_1475f <= 1000
 SUCCt2_2__TaxaB <= 1000
 SUCCt2_2__TaxaB_reverse_4bdaf = 0
 SUCCt3__TaxaB <= 1000
 SUCCt3__TaxaB_reverse_f38cd = 0
 SUCDi__TaxaB <= 1000
 SUCDi__TaxaB_reverse_c5ffa = 0
 SUCOAS__TaxaB <= 1000
 SUCOAS__TaxaB_reverse_3f8e7 <= 1000
 TALA__TaxaB <= 1000
 TALA__TaxaB_reverse_050f7 <= 1000
 THD2__TaxaB <= 1000
 THD2__TaxaB_reverse_c7723 = 0
 TKT1__TaxaB <= 1000
 TKT1__TaxaB_reverse_d1c70 <= 1000
 TKT2__TaxaB <= 1000
 TKT2__TaxaB_reverse_39177 <= 1000
 TPI__TaxaB <= 1000
 TPI__TaxaB_reverse_b9846 <= 1000
 ACALD__TaxaC <= 1000
 ACALD__TaxaC_reverse_b9484 <= 1000
 ACALDt__TaxaC <= 1000
 ACALDt__TaxaC_reverse_d1aa6 <= 1000
 ACKr__TaxaC <= 1000
 ACKr__TaxaC_reverse_be2cc <= 1000
 ACONTa__TaxaC <= 1000
 ACONTa__TaxaC_reverse_21b27 <= 1000
 ACONTb__TaxaC <= 1000
 ACONTb__TaxaC_reverse_9f088 <= 1000
 ACt2r__TaxaC <= 1000
 ACt2r__TaxaC_reverse_b30b4 <= 1000
 ADK1__TaxaC <= 1000
 ADK1__TaxaC_reverse_0a628 <= 1000
 AKGDH__TaxaC <= 1000
 AKGDH__TaxaC_reverse_a3c9d = 0
 AKGt2r__TaxaC <= 1000
 AKGt2r__TaxaC_reverse_72a4e <= 1000
 ALCD2x__TaxaC <= 1000
 ALCD2x__TaxaC_reverse_63b1c <= 1000
 8.39 <= ATPM__TaxaC <= 1000
 ATPM__TaxaC_reverse_d3438 = 0
 ATPS4r__TaxaC <= 1000
 ATPS4r__TaxaC_reverse_edfcc <= 1000
 biomass_ecoli__TaxaC <= 1000
 biomass_ecoli__TaxaC_reverse_41d18 = 0
 CO2t__TaxaC <= 1000
 CO2t__TaxaC_reverse_3ea53 <= 1000
 CS__TaxaC <= 1000
 CS__TaxaC_reverse_162ab = 0
 CYTBD__TaxaC <= 1000
 CYTBD__TaxaC_reverse_b5e26 = 0
 D_LACt2__TaxaC <= 1000
 D_LACt2__TaxaC_reverse_fb45c <= 1000
 ENO__TaxaC <= 1000
 ENO__TaxaC_reverse_e5c24 <= 1000
 ETOHt2r__TaxaC <= 1000
 ETOHt2r__TaxaC_reverse_5a87c <= 1000
 IEX_ac[u]__TaxaC <= 1000
 IEX_ac[u]__TaxaC_reverse_82616 <= 1000
 IEX_acald[u]__TaxaC <= 1000
 IEX_acald[u]__TaxaC_reverse_79433 <= 1000
 IEX_akg[u]__TaxaC <= 1000
 IEX_akg[u]__TaxaC_reverse_2d021 <= 1000
 IEX_co2[u]__TaxaC <= 1000
 IEX_co2[u]__TaxaC_reverse_f8836 <= 1000
 EX[u]toh[u]__TaxaC <= 1000
 EX[u]toh[u]__TaxaC_reverse_65089 <= 1000
 IEX_for[u]__TaxaC <= 1000
 IEX_for[u]__TaxaC_reverse_89537 <= 1000
 IEX_fru[u]__TaxaC <= 1000
 IEX_fru[u]__TaxaC_reverse_d1617 <= 1000
 IEX_fum[u]__TaxaC <= 1000
 IEX_fum[u]__TaxaC_reverse_66a1c <= 1000
 IEX_glc__D[u]__TaxaC <= 1000
 IEX_glc__D[u]__TaxaC_reverse_4dbc9 <= 1000
 IEX_gln__L[u]__TaxaC <= 1000
 IEX_gln__L[u]__TaxaC_reverse_d4f65 <= 1000
 IEX_glu__L[u]__TaxaC <= 1000
 IEX_glu__L[u]__TaxaC_reverse_b053a <= 1000
 IEX_h[u]__TaxaC <= 1000
 IEX_h[u]__TaxaC_reverse_b39fd <= 1000
 IEX_h2o[u]__TaxaC <= 1000
 IEX_h2o[u]__TaxaC_reverse_f3c2c <= 1000
 IEX_lac__D[u]__TaxaC <= 1000
 IEX_lac__D[u]__TaxaC_reverse_baf61 <= 1000
 IEX_mal__L[u]__TaxaC <= 1000
 IEX_mal__L[u]__TaxaC_reverse_2b6d1 <= 1000
 IEX_nh4[u]__TaxaC <= 1000
 IEX_nh4[u]__TaxaC_reverse_486f8 <= 1000
 IEX_o2[u]__TaxaC <= 1000
 IEX_o2[u]__TaxaC_reverse_9d59a <= 1000
 IEX_pi[u]__TaxaC <= 1000
 IEX_pi[u]__TaxaC_reverse_acfe9 <= 1000
 IEX_pyr[u]__TaxaC <= 1000
 IEX_pyr[u]__TaxaC_reverse_74aed <= 1000
 IEX_succ[u]__TaxaC <= 1000
 IEX_succ[u]__TaxaC_reverse_0ba26 <= 1000
 FBA__TaxaC <= 1000
 FBA__TaxaC_reverse_50114 <= 1000
 FBP__TaxaC <= 1000
 FBP__TaxaC_reverse_27fed = 0
 FORt2__TaxaC <= 1000
 FORt2__TaxaC_reverse_c2407 = 0
 FORti__TaxaC <= 1000
 FORti__TaxaC_reverse_ef905 = 0
 FRD7__TaxaC <= 1000
 FRD7__TaxaC_reverse_eec45 = 0
 FRUpts2__TaxaC <= 1000
 FRUpts2__TaxaC_reverse_01b00 = 0
 FUM__TaxaC <= 1000
 FUM__TaxaC_reverse_0d074 <= 1000
 FUMt2_2__TaxaC <= 1000
 FUMt2_2__TaxaC_reverse_4db1b = 0
 G6PDH2r__TaxaC <= 1000
 G6PDH2r__TaxaC_reverse_ab3e1 <= 1000
 GAPD__TaxaC <= 1000
 GAPD__TaxaC_reverse_17a1a <= 1000
 GLCpts__TaxaC <= 1000
 GLCpts__TaxaC_reverse_d4b50 = 0
 GLNS__TaxaC <= 1000
 GLNS__TaxaC_reverse_12565 = 0
 GLNabc__TaxaC <= 1000
 GLNabc__TaxaC_reverse_904d0 = 0
 GLUDy__TaxaC <= 1000
 GLUDy__TaxaC_reverse_990ca <= 1000
 GLUN__TaxaC <= 1000
 GLUN__TaxaC_reverse_d2cdc = 0
 GLUSy__TaxaC <= 1000
 GLUSy__TaxaC_reverse_84635 = 0
 GLUt2r__TaxaC <= 1000
 GLUt2r__TaxaC_reverse_b3974 <= 1000
 GND__TaxaC <= 1000
 GND__TaxaC_reverse_87f3c = 0
 H2Ot__TaxaC <= 1000
 H2Ot__TaxaC_reverse_3f623 <= 1000
 ICDHyr__TaxaC <= 1000
 ICDHyr__TaxaC_reverse_5612a <= 1000
 ICL__TaxaC <= 1000
 ICL__TaxaC_reverse_6b407 = 0
 LDH_D__TaxaC <= 1000
 LDH_D__TaxaC_reverse_56e41 <= 1000
 MALS__TaxaC <= 1000
 MALS__TaxaC_reverse_da103 = 0
 MALt2_2__TaxaC <= 1000
 MALt2_2__TaxaC_reverse_680c0 = 0
 MDH__TaxaC <= 1000
 MDH__TaxaC_reverse_5bd24 <= 1000
 ME1__TaxaC <= 1000
 ME1__TaxaC_reverse_b9203 = 0
 ME2__TaxaC <= 1000
 ME2__TaxaC_reverse_cf032 = 0
 NADH16__TaxaC <= 1000
 NADH16__TaxaC_reverse_fb2ba = 0
 NADTRHD__TaxaC <= 1000
 NADTRHD__TaxaC_reverse_ac758 = 0
 NH4t__TaxaC <= 1000
 NH4t__TaxaC_reverse_b46fe <= 1000
 O2t__TaxaC <= 1000
 O2t__TaxaC_reverse_6c905 <= 1000
 PDH__TaxaC <= 1000
 PDH__TaxaC_reverse_87edc = 0
 PFK__TaxaC <= 1000
 PFK__TaxaC_reverse_a987d = 0
 PFL__TaxaC <= 1000
 PFL__TaxaC_reverse_3769b = 0
 PGI__TaxaC <= 1000
 PGI__TaxaC_reverse_b32d8 <= 1000
 PGK__TaxaC <= 1000
 PGK__TaxaC_reverse_aaf6a <= 1000
 PGL__TaxaC <= 1000
 PGL__TaxaC_reverse_0181e = 0
 PGM__TaxaC <= 1000
 PGM__TaxaC_reverse_cfa90 <= 1000
 PIt2r__TaxaC <= 1000
 PIt2r__TaxaC_reverse_206b5 <= 1000
 PPC__TaxaC <= 1000
 PPC__TaxaC_reverse_b4752 = 0
 PPCK__TaxaC <= 1000
 PPCK__TaxaC_reverse_400f8 = 0
 PPS__TaxaC <= 1000
 PPS__TaxaC_reverse_dba7f = 0
 PTAr__TaxaC <= 1000
 PTAr__TaxaC_reverse_48ac3 <= 1000
 PYK__TaxaC <= 1000
 PYK__TaxaC_reverse_95865 = 0
 PYRt2__TaxaC <= 1000
 PYRt2__TaxaC_reverse_aadfb <= 1000
 RPE__TaxaC <= 1000
 RPE__TaxaC_reverse_8923a <= 1000
 RPI__TaxaC <= 1000
 RPI__TaxaC_reverse_2ab80 <= 1000
 SUCCt2_2__TaxaC <= 1000
 SUCCt2_2__TaxaC_reverse_572ac = 0
 SUCCt3__TaxaC <= 1000
 SUCCt3__TaxaC_reverse_a206f = 0
 SUCDi__TaxaC <= 1000
 SUCDi__TaxaC_reverse_10c32 = 0
 SUCOAS__TaxaC <= 1000
 SUCOAS__TaxaC_reverse_e07eb <= 1000
 TALA__TaxaC <= 1000
 TALA__TaxaC_reverse_8e9b1 <= 1000
 THD2__TaxaC <= 1000
 THD2__TaxaC_reverse_a7f74 = 0
 TKT1__TaxaC <= 1000
 TKT1__TaxaC_reverse_61d41 <= 1000
 TKT2__TaxaC <= 1000
 TKT2__TaxaC_reverse_f5ba9 <= 1000
 TPI__TaxaC <= 1000
 TPI__TaxaC_reverse_6b386 <= 1000
 ACALD__TaxaD <= 1000
 ACALD__TaxaD_reverse_6d0be <= 1000
 ACALDt__TaxaD <= 1000
 ACALDt__TaxaD_reverse_0a331 <= 1000
 ACKr__TaxaD <= 1000
 ACKr__TaxaD_reverse_f63fc <= 1000
 ACONTa__TaxaD <= 1000
 ACONTa__TaxaD_reverse_013c6 <= 1000
 ACONTb__TaxaD <= 1000
 ACONTb__TaxaD_reverse_799dc <= 1000
 ACt2r__TaxaD <= 1000
 ACt2r__TaxaD_reverse_36c26 <= 1000
 ADK1__TaxaD <= 1000
 ADK1__TaxaD_reverse_af0f9 <= 1000
 AKGDH__TaxaD <= 1000
 AKGDH__TaxaD_reverse_049a5 = 0
 AKGt2r__TaxaD <= 1000
 AKGt2r__TaxaD_reverse_2bddb <= 1000
 ALCD2x__TaxaD <= 1000
 ALCD2x__TaxaD_reverse_f3fcd <= 1000
 8.39 <= ATPM__TaxaD <= 1000
 ATPM__TaxaD_reverse_39785 = 0
 ATPS4r__TaxaD <= 1000
 ATPS4r__TaxaD_reverse_5028a <= 1000
 biomass_ecoli__TaxaD <= 1000
 biomass_ecoli__TaxaD_reverse_71a87 = 0
 CO2t__TaxaD <= 1000
 CO2t__TaxaD_reverse_f34ba <= 1000
 CS__TaxaD <= 1000
 CS__TaxaD_reverse_aaf14 = 0
 CYTBD__TaxaD <= 1000
 CYTBD__TaxaD_reverse_4bb10 = 0
 D_LACt2__TaxaD <= 1000
 D_LACt2__TaxaD_reverse_b0332 <= 1000
 ENO__TaxaD <= 1000
 ENO__TaxaD_reverse_7d37f <= 1000
 ETOHt2r__TaxaD <= 1000
 ETOHt2r__TaxaD_reverse_9310d <= 1000
 IEX_ac[u]__TaxaD <= 1000
 IEX_ac[u]__TaxaD_reverse_63ca5 <= 1000
 IEX_acald[u]__TaxaD <= 1000
 IEX_acald[u]__TaxaD_reverse_1bc86 <= 1000
 IEX_akg[u]__TaxaD <= 1000
 IEX_akg[u]__TaxaD_reverse_a37cb <= 1000
 IEX_co2[u]__TaxaD <= 1000
 IEX_co2[u]__TaxaD_reverse_ce75e <= 1000
 EX[u]toh[u]__TaxaD <= 1000
 EX[u]toh[u]__TaxaD_reverse_1d69d <= 1000
 IEX_for[u]__TaxaD <= 1000
 IEX_for[u]__TaxaD_reverse_0a866 <= 1000
 IEX_fru[u]__TaxaD <= 1000
 IEX_fru[u]__TaxaD_reverse_8c9bf <= 1000
 IEX_fum[u]__TaxaD <= 1000
 IEX_fum[u]__TaxaD_reverse_3415d <= 1000
 IEX_glc__D[u]__TaxaD <= 1000
 IEX_glc__D[u]__TaxaD_reverse_405dc <= 1000
 IEX_gln__L[u]__TaxaD <= 1000
 IEX_gln__L[u]__TaxaD_reverse_23702 <= 1000
 IEX_glu__L[u]__TaxaD <= 1000
 IEX_glu__L[u]__TaxaD_reverse_db289 <= 1000
 IEX_h[u]__TaxaD <= 1000
 IEX_h[u]__TaxaD_reverse_3e6c7 <= 1000
 IEX_h2o[u]__TaxaD <= 1000
 IEX_h2o[u]__TaxaD_reverse_c1c75 <= 1000
 IEX_lac__D[u]__TaxaD <= 1000
 IEX_lac__D[u]__TaxaD_reverse_115c3 <= 1000
 IEX_mal__L[u]__TaxaD <= 1000
 IEX_mal__L[u]__TaxaD_reverse_36f62 <= 1000
 IEX_nh4[u]__TaxaD <= 1000
 IEX_nh4[u]__TaxaD_reverse_a85a3 <= 1000
 IEX_o2[u]__TaxaD <= 1000
 IEX_o2[u]__TaxaD_reverse_a1b64 <= 1000
 IEX_pi[u]__TaxaD <= 1000
 IEX_pi[u]__TaxaD_reverse_925d1 <= 1000
 IEX_pyr[u]__TaxaD <= 1000
 IEX_pyr[u]__TaxaD_reverse_0e464 <= 1000
 IEX_succ[u]__TaxaD <= 1000
 IEX_succ[u]__TaxaD_reverse_34a9d <= 1000
 FBA__TaxaD <= 1000
 FBA__TaxaD_reverse_4068a <= 1000
 FBP__TaxaD <= 1000
 FBP__TaxaD_reverse_5bb38 = 0
 FORt2__TaxaD <= 1000
 FORt2__TaxaD_reverse_8e3d1 = 0
 FORti__TaxaD <= 1000
 FORti__TaxaD_reverse_a472d = 0
 FRD7__TaxaD <= 1000
 FRD7__TaxaD_reverse_2b38f = 0
 FRUpts2__TaxaD <= 1000
 FRUpts2__TaxaD_reverse_61ab1 = 0
 FUM__TaxaD <= 1000
 FUM__TaxaD_reverse_3a7c0 <= 1000
 FUMt2_2__TaxaD <= 1000
 FUMt2_2__TaxaD_reverse_d9367 = 0
 G6PDH2r__TaxaD <= 1000
 G6PDH2r__TaxaD_reverse_a3f2e <= 1000
 GAPD__TaxaD <= 1000
 GAPD__TaxaD_reverse_20540 <= 1000
 GLCpts__TaxaD <= 1000
 GLCpts__TaxaD_reverse_e25bf = 0
 GLNS__TaxaD <= 1000
 GLNS__TaxaD_reverse_5ed39 = 0
 GLNabc__TaxaD <= 1000
 GLNabc__TaxaD_reverse_85438 = 0
 GLUDy__TaxaD <= 1000
 GLUDy__TaxaD_reverse_203e2 <= 1000
 GLUN__TaxaD <= 1000
 GLUN__TaxaD_reverse_4524e = 0
 GLUSy__TaxaD <= 1000
 GLUSy__TaxaD_reverse_f2768 = 0
 GLUt2r__TaxaD <= 1000
 GLUt2r__TaxaD_reverse_cedc6 <= 1000
 GND__TaxaD <= 1000
 GND__TaxaD_reverse_efd8f = 0
 H2Ot__TaxaD <= 1000
 H2Ot__TaxaD_reverse_37165 <= 1000
 ICDHyr__TaxaD <= 1000
 ICDHyr__TaxaD_reverse_f264a <= 1000
 ICL__TaxaD <= 1000
 ICL__TaxaD_reverse_715d6 = 0
 LDH_D__TaxaD <= 1000
 LDH_D__TaxaD_reverse_49604 <= 1000
 MALS__TaxaD <= 1000
 MALS__TaxaD_reverse_e8b65 = 0
 MALt2_2__TaxaD <= 1000
 MALt2_2__TaxaD_reverse_cf2d6 = 0
 MDH__TaxaD <= 1000
 MDH__TaxaD_reverse_b8d94 <= 1000
 ME1__TaxaD <= 1000
 ME1__TaxaD_reverse_e55a1 = 0
 ME2__TaxaD <= 1000
 ME2__TaxaD_reverse_d911c = 0
 NADH16__TaxaD <= 1000
 NADH16__TaxaD_reverse_a5f47 = 0
 NADTRHD__TaxaD <= 1000
 NADTRHD__TaxaD_reverse_22e5b = 0
 NH4t__TaxaD <= 1000
 NH4t__TaxaD_reverse_27ea9 <= 1000
 O2t__TaxaD <= 1000
 O2t__TaxaD_reverse_169fb <= 1000
 PDH__TaxaD <= 1000
 PDH__TaxaD_reverse_a3fa6 = 0
 PFK__TaxaD <= 1000
 PFK__TaxaD_reverse_846b7 = 0
 PFL__TaxaD <= 1000
 PFL__TaxaD_reverse_a5dcf = 0
 PGI__TaxaD <= 1000
 PGI__TaxaD_reverse_6ff0e <= 1000
 PGK__TaxaD <= 1000
 PGK__TaxaD_reverse_15b8f <= 1000
 PGL__TaxaD <= 1000
 PGL__TaxaD_reverse_91fdc = 0
 PGM__TaxaD <= 1000
 PGM__TaxaD_reverse_732f2 <= 1000
 PIt2r__TaxaD <= 1000
 PIt2r__TaxaD_reverse_93097 <= 1000
 PPC__TaxaD <= 1000
 PPC__TaxaD_reverse_b4029 = 0
 PPCK__TaxaD <= 1000
 PPCK__TaxaD_reverse_30ef8 = 0
 PPS__TaxaD <= 1000
 PPS__TaxaD_reverse_fdd65 = 0
 PTAr__TaxaD <= 1000
 PTAr__TaxaD_reverse_bb1d7 <= 1000
 PYK__TaxaD <= 1000
 PYK__TaxaD_reverse_70730 = 0
 PYRt2__TaxaD <= 1000
 PYRt2__TaxaD_reverse_20065 <= 1000
 RPE__TaxaD <= 1000
 RPE__TaxaD_reverse_82d61 <= 1000
 RPI__TaxaD <= 1000
 RPI__TaxaD_reverse_d146e <= 1000
 SUCCt2_2__TaxaD <= 1000
 SUCCt2_2__TaxaD_reverse_29f5f = 0
 SUCCt3__TaxaD <= 1000
 SUCCt3__TaxaD_reverse_282b4 = 0
 SUCDi__TaxaD <= 1000
 SUCDi__TaxaD_reverse_35a39 = 0
 SUCOAS__TaxaD <= 1000
 SUCOAS__TaxaD_reverse_21a9f <= 1000
 TALA__TaxaD <= 1000
 TALA__TaxaD_reverse_d504b <= 1000
 THD2__TaxaD <= 1000
 THD2__TaxaD_reverse_7e9dc = 0
 TKT1__TaxaD <= 1000
 TKT1__TaxaD_reverse_f64cc <= 1000
 TKT2__TaxaD <= 1000
 TKT2__TaxaD_reverse_68f4a <= 1000
 TPI__TaxaD <= 1000
 TPI__TaxaD_reverse_1efde <= 1000
 EX_microbeBiomass[fe] <= 1e+06
 EX_microbeBiomass[fe]_reverse_784f8 <= 1000
 UFEt_microbeBiomass <= 1e+06
 UFEt_microbeBiomass_reverse_e9f5a = 0
 0.4 <= communityBiomass <= 1
 communityBiomass_reverse_5ad9d = 0
End
