$ root ../data/ORCA_MC_summary_all_10Apr2018.root
summary->Draw("fMC_energy>>h1(100,0,100)", "fRDF_track_score >= 0.9")
c1->SaveAs("th1d_mc_energy_q_lt_0.9.pdf")
