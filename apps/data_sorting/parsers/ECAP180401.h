//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 22 13:20:36 2019 by ROOT version 6.10/02
// from TTree PID/
// found on file: irodsdata/pid_output_atm_neutrino_atm_muon_pure_noise_shiftedVertexEventSelection_180401.root
//////////////////////////////////////////////////////////

#ifndef ECAP180401_h
#define ECAP180401_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ECAP180401 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        distance_dusj_recolns;
   Double_t        recolns_pos_r;
   Double_t        dusj_theta;
   Double_t        recolns_dir_xyz;
   Double_t        recolns_shifted30m_pos_z;
   Double_t        dusj_log10_energy_corrected;
   Double_t        recolns_n_hits_used;
   Double_t        gandalf_phi;
   Double_t        gandalf_is_selected_without_containment;
   Double_t        dusj_sumMeasuredPMThits;
   Double_t        dusj_MuonSuppression_enoughHits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_energy;
   Double_t        recolns_shifted30m_pos_y;
   Double_t        recolns_quality;
   Double_t        _crossValidationTrainingMask_noise_score;
   Double_t        distance_truth_gandalf;
   Double_t        dusj_L0AroundL1HitSelection__weight_withoutSelfSquared;
   Double_t        dusj_shifted30m_pos_x;
   Double_t        dusj_passes_anti_noise_cut;
   Double_t        dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80;
   Double_t        dusj_phi;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge;
   Double_t        gandalf_spread_pos_x_median;
   Double_t        mc_id;
   Double_t        dusj_L0AroundL1HitSelection_weight_withoutSelfSquared;
   Double_t        gandalf_spread_pos_y_iqr;
   Double_t        dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference;
   Double_t        dusj_llhBestSinglePMTperDOM__sum_forNoSignal;
   Double_t        dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_y;
   Double_t        dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_x;
   Double_t        dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_z;
   Double_t        log10_energy;
   Double_t        gandalf_spread_dir_z_mad;
   Double_t        gandalf_is_selected;
   Double_t        dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N;
   Double_t        noise_score;
   Double_t        gandalf_jstart_length;
   Double_t        dusj_Nom;
   Double_t        cosangle_gandalf_recolns;
   Double_t        gandalf_cosangle_dir_pos;
   Double_t        gandalf_lambda;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_phi;
   Double_t        selected_in_ECAP_PID;
   Double_t        dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber;
   Double_t        dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference;
   Double_t        gandalf_containment;
   Double_t        gandalf_spread_dir_x_mean;
   Double_t        pos_z;
   Double_t        pos_x;
   Double_t        pos_y;
   Double_t        gandalf_beta1;
   Double_t        dusj_energy;
   Double_t        gandalf_beta0;
   Double_t        pos_r;
   Double_t        dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75;
   Double_t        gandalf_spread_dir_y_std;
   Double_t        recolns_loose_is_selected;
   Double_t        recolns_cosangle_dir_pos;
   Double_t        dusj_Fork_muonSuppression_decision;
   Double_t        gandalf_spread_pos_z_std;
   Double_t        gandalf_spread_pos_z_median;
   Double_t        cosangle_dusj_gandalf;
   Double_t        gandalf_jenergy_energy;
   Double_t        livetime_sec;
   Double_t        recolns_is_good;
   Double_t        cos_theta;
   Double_t        is_cc;
   Double_t        gandalf_n_hits;
   Double_t        n_files_gen;
   Double_t        gandalf_spread_beta0_iqr;
   Double_t        energy;
   Double_t        dusj_Npmt;
   Double_t        Emuon_from_neutrino_interaction;
   Double_t        recolns_shifted30m_pos_x;
   Double_t        recolns_log10_energy_muon;
   Double_t        is_neutrino;
   Double_t        gandalf_pos_y;
   Double_t        dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;
   Double_t        n_events_gen;
   Double_t        dir_xyz;
   Double_t        gandalf_spread_n_hits_mean;
   Double_t        dusj_DifferencesFirstAndSecondVertexFit_deltaTime;
   Double_t        recolns_sigma2_theta;
   Double_t        recolns_phi;
   Double_t        trigger_mask;
   Double_t        gandalf_energy_corrected;
   Double_t        dir_x;
   Double_t        recolns_pos_y;
   Double_t        cosangle_truth_recolns;
   Double_t        recolns_pos_z;
   Double_t        gandalf_spread_pos_y_median;
   Double_t        gandalf_spread_dir_x_iqr;
   Double_t        dusj_llhSinglePMT_sum;
   Double_t        _crossValidationTrainingMask_muon_score;
   Double_t        gandalf_spread_dir_z_iqr;
   Double_t        dusj_is_selected_without_anti_noise_cut;
   Double_t        utc_nanoseconds;
   Double_t        gandalf_spread_beta0_median;
   Double_t        gandalf_pos_r;
   Double_t        recolns_pos_xyz;
   Double_t        dusj_llh_sum;
   Double_t        bjorkeny;
   Double_t        theta;
   Double_t        dusj_llhSinglePMT_sum_forNoSignal;
   Double_t        gandalf_spread_chi2_std;
   Double_t        interaction_channel;
   Double_t        gandalf_pos_z;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge;
   Double_t        gandalf_pos_x;
   Double_t        dusj_targetEnergy;
   Double_t        dir_y;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean;
   Double_t        dusj_energyErrorDown_bestLLH;
   Double_t        dir_z;
   Double_t        dusj_EventID;
   Double_t        gandalf_spread_pos_z_mean;
   Double_t        dir_r;
   Double_t        gandalf_spread_n_hits_std;
   Double_t        dusj_best_SecondDusjOrcaVertexFit_FitResult_time;
   Double_t        dusj_shifted30m_pos_y;
   Double_t        dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber;
   Double_t        dusj_shifted30m_pos_z;
   Double_t        pos_xyz;
   Double_t        dusj_zenith;
   Double_t        gandalf_is_good;
   Double_t        dusj_Npmt_maxDeltaT10ns;
   Double_t        gandalf_shifted30m_pos_x;
   Double_t        dusj_dir_xyz;
   Double_t        gandalf_spread_chi2_median;
   Double_t        gandalf_dir_y;
   Double_t        gandalf_dir_x;
   Double_t        gandalf_log10_energy_corrected;
   Double_t        gandalf_dir_z;
   Double_t        distancePointTrack_dusj_truth;
   Double_t        gandalf_loose_is_selected;
   Double_t        gandalf_dir_r;
   Double_t        gandalf_spread_dir_z_median;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean;
   Double_t        length;
   Double_t        gandalf_spread_dir_z_mean;
   Double_t        gandalf_spread_beta1_mean;
   Double_t        gandalf_dir_xyz;
   Double_t        dusj_pos_xyz;
   Double_t        gandalf_spread_pos_y_mean;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct;
   Double_t        recolns_is_selected_without_containment;
   Double_t        dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N;
   Double_t        gandalf_spread_pos_x_mean;
   Double_t        gandalf_cos_theta;
   Double_t        distance_truth_recolns;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_time;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient;
   Double_t        gandalf_spread_pos_z_mad;
   Double_t        gandalf_zenith;
   Double_t        run_id;
   Double_t        dusj_energy_corrected;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY;
   Double_t        dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_theta;
   Double_t        tracklike_from_Emuon;
   Double_t        dusj_is_selected;
   Double_t        gandalf_spread_chi2_mean;
   Double_t        dusj_DifferencesFirstAndSecondVertexFit_distance;
   Double_t        dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;
   Double_t        dusj_best_FirstDusjOrcaVertexFit_FitResult_time;
   Double_t        dusj_azimuth;
   Double_t        recolns_theta;
   Double_t        recolns_zenith;
   Double_t        gandalf_time;
   Double_t        gandalf_spread_dir_x_std;
   Double_t        gandalf_containment_shifted;
   Double_t        recolns_is_selected;
   Double_t        gandalf_spread_beta1_iqr;
   Double_t        dusj_cosangle_dir_pos;
   Double_t        gandalf_spread_dir_x_median;
   Double_t        dusj_MuonSuppression_decision;
   Double_t        recolns_sigma2_phi;
   Double_t        gandalf_jstart_npe_mip_total;
   Double_t        type;
   Double_t        dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;
   Double_t        dusj_L0AroundL1HitSelection_weight;
   Double_t        recolns_dir_x;
   Double_t        recolns_dir_y;
   Double_t        log_distance4Dshortest_dusj_gandalf;
   Double_t        recolns_dir_r;
   Double_t        cosangle_dusj_recolns;
   Double_t        dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber;
   Double_t        gandalf_spread_beta0_std;
   Double_t        gandalf_spread_n_hits_median;
   Double_t        trigger_counter;
   Double_t        recolns_log10_energy_neutrino;
   Double_t        distancePointTrack_gandalf_truth;
   Double_t        dusj_time;
   Double_t        gandalf_spread_pos_y_mad;
   Double_t        gandalf_chi2;
   Double_t        cos_zenith;
   Double_t        gandalf_jstart_npe_mip;
   Double_t        gandalf_azimuth;
   Double_t        gandalf_spread_dir_z_std;
   Double_t        recolns_energy_neutrino;
   Double_t        dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_z;
   Double_t        track_score;
   Double_t        dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_x;
   Double_t        dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_y;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_z;
   Double_t        mc_t;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_x;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_y;
   Double_t        muon_score;
   Double_t        dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;
   Double_t        gandalf_theta;
   Double_t        dusj_dir_r;
   Double_t        gandalf_n_iter;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean;
   Double_t        gandalf_spread_n_hits_iqr;
   Double_t        recolns_energy_muon;
   Double_t        recolns_containment;
   Double_t        weight_one_year;
   Double_t        distance_dusj_gandalf;
   Double_t        det_id;
   Double_t        recolns_dir_z;
   Double_t        cosangle_truth_dusj;
   Double_t        distance_gandalf_recolns;
   Double_t        zenith;
   Double_t        gandalf_pos_xyz;
   Double_t        dusj_sumMeasuredOMhits;
   Double_t        utc_seconds;
   Double_t        gandalf_spread_pos_y_std;
   Double_t        gandalf_spread_pos_x_std;
   Double_t        dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N;
   Double_t        dusj_sumExpOMhits;
   Double_t        Erange_min;
   Double_t        dusj_pos_r;
   Double_t        gandalf_cos_zenith;
   Double_t        gandalf_spread_beta1_median;
   Double_t        overlays;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits;
   Double_t        id;
   Double_t        dusj_pos_x;
   Double_t        distance_truth_dusj;
   Double_t        dusj_energyErrorUp_bestLLH;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber;
   Double_t        event_id;
   Double_t        dusj_pos_z;
   Double_t        cosangle_truth_gandalf;
   Double_t        dusj_is_good;
   Double_t        gandalf_spread_pos_x_mad;
   Double_t        Erange_max;
   Double_t        gandalf_spread_pos_x_iqr;
   Double_t        gandalf_spread_dir_y_iqr;
   Double_t        gandalf_spread_beta1_std;
   Double_t        dusj_cos_zenith;
   Double_t        gandalf_spread_beta0_mean;
   Double_t        recolns_pos_x;
   Double_t        dusj_pos_y;
   Double_t        dusj_energy_bestLLH;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean;
   Double_t        phi;
   Double_t        dusj_llhBestSinglePMTperDOM_sum;
   Double_t        weight_w3;
   Double_t        weight_w2;
   Double_t        weight_w1;
   Double_t        dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct;
   Double_t        azimuth;
   Double_t        recolns_sin_theta;
   Double_t        gandalf_spread_beta1_mad;
   Double_t        gandalf_spread_chi2_mad;
   Double_t        gandalf_jenergy_chi2;
   Double_t        gandalf_spread_chi2_iqr;
   Double_t        recolns_containment_shifted;
   Double_t        recolns_bjorken_y;
   Double_t        recolns_cos_theta;
   Double_t        recolns_beta;
   Double_t        gandalf_spread_pos_z_iqr;
   Double_t        dusj_cos_theta;
   Double_t        dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber;
   Double_t        gandalf_shifted30m_pos_y;
   Double_t        frame_index;
   Double_t        gandalf_shifted30m_pos_z;
   Double_t        dusj_sumExpFromPoissonOMhits;
   Double_t        gandalf_spread_dir_x_mad;
   Double_t        dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75;
   Double_t        _crossValidationTrainingMask_track_score;
   Double_t        gandalf_spread_dir_y_mean;
   Double_t        gandalf_spread_beta0_mad;
   Double_t        dusj_multiplicity;
   Double_t        gandalf_spread_n_hits_mad;
   Double_t        log_distance_dusj_gandalfForDusjTime;
   Double_t        dusj_MuonSuppression_deltaTresQ20Q80;
   Double_t        gandalf_spread_dir_y_median;
   Double_t        time;
   Double_t        recolns_azimuth;
   Double_t        dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference;
   Double_t        gandalf_spread_dir_y_mad;
   Double_t        recolns_cos_zenith;
   Double_t        dusj_dir_y;
   Double_t        dusj_dir_x;
   Double_t        dusj_dir_z;

   // List of branches
   TBranch        *b_distance_dusj_recolns;   //!
   TBranch        *b_recolns_pos_r;   //!
   TBranch        *b_dusj_theta;   //!
   TBranch        *b_recolns_dir_xyz;   //!
   TBranch        *b_recolns_shifted30m_pos_z;   //!
   TBranch        *b_dusj_log10_energy_corrected;   //!
   TBranch        *b_recolns_n_hits_used;   //!
   TBranch        *b_gandalf_phi;   //!
   TBranch        *b_gandalf_is_selected_without_containment;   //!
   TBranch        *b_dusj_sumMeasuredPMThits;   //!
   TBranch        *b_dusj_MuonSuppression_enoughHits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_energy;   //!
   TBranch        *b_recolns_shifted30m_pos_y;   //!
   TBranch        *b_recolns_quality;   //!
   TBranch        *b__crossValidationTrainingMask_noise_score;   //!
   TBranch        *b_distance_truth_gandalf;   //!
   TBranch        *b_dusj_L0AroundL1HitSelection__weight_withoutSelfSquared;   //!
   TBranch        *b_dusj_shifted30m_pos_x;   //!
   TBranch        *b_dusj_passes_anti_noise_cut;   //!
   TBranch        *b_dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80;   //!
   TBranch        *b_dusj_phi;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge;   //!
   TBranch        *b_gandalf_spread_pos_x_median;   //!
   TBranch        *b_mc_id;   //!
   TBranch        *b_dusj_L0AroundL1HitSelection_weight_withoutSelfSquared;   //!
   TBranch        *b_gandalf_spread_pos_y_iqr;   //!
   TBranch        *b_dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference;   //!
   TBranch        *b_dusj_llhBestSinglePMTperDOM__sum_forNoSignal;   //!
   TBranch        *b_dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_y;   //!
   TBranch        *b_dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_x;   //!
   TBranch        *b_dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_z;   //!
   TBranch        *b_log10_energy;   //!
   TBranch        *b_gandalf_spread_dir_z_mad;   //!
   TBranch        *b_gandalf_is_selected;   //!
   TBranch        *b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N;   //!
   TBranch        *b_noise_score;   //!
   TBranch        *b_gandalf_jstart_length;   //!
   TBranch        *b_dusj_Nom;   //!
   TBranch        *b_cosangle_gandalf_recolns;   //!
   TBranch        *b_gandalf_cosangle_dir_pos;   //!
   TBranch        *b_gandalf_lambda;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_phi;   //!
   TBranch        *b_selected_in_ECAP_PID;   //!
   TBranch        *b_dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber;   //!
   TBranch        *b_dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference;   //!
   TBranch        *b_gandalf_containment;   //!
   TBranch        *b_gandalf_spread_dir_x_mean;   //!
   TBranch        *b_pos_z;   //!
   TBranch        *b_pos_x;   //!
   TBranch        *b_pos_y;   //!
   TBranch        *b_gandalf_beta1;   //!
   TBranch        *b_dusj_energy;   //!
   TBranch        *b_gandalf_beta0;   //!
   TBranch        *b_pos_r;   //!
   TBranch        *b_dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75;   //!
   TBranch        *b_gandalf_spread_dir_y_std;   //!
   TBranch        *b_recolns_loose_is_selected;   //!
   TBranch        *b_recolns_cosangle_dir_pos;   //!
   TBranch        *b_dusj_Fork_muonSuppression_decision;   //!
   TBranch        *b_gandalf_spread_pos_z_std;   //!
   TBranch        *b_gandalf_spread_pos_z_median;   //!
   TBranch        *b_cosangle_dusj_gandalf;   //!
   TBranch        *b_gandalf_jenergy_energy;   //!
   TBranch        *b_livetime_sec;   //!
   TBranch        *b_recolns_is_good;   //!
   TBranch        *b_cos_theta;   //!
   TBranch        *b_is_cc;   //!
   TBranch        *b_gandalf_n_hits;   //!
   TBranch        *b_n_files_gen;   //!
   TBranch        *b_gandalf_spread_beta0_iqr;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_dusj_Npmt;   //!
   TBranch        *b_Emuon_from_neutrino_interaction;   //!
   TBranch        *b_recolns_shifted30m_pos_x;   //!
   TBranch        *b_recolns_log10_energy_muon;   //!
   TBranch        *b_is_neutrino;   //!
   TBranch        *b_gandalf_pos_y;   //!
   TBranch        *b_dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;   //!
   TBranch        *b_n_events_gen;   //!
   TBranch        *b_dir_xyz;   //!
   TBranch        *b_gandalf_spread_n_hits_mean;   //!
   TBranch        *b_dusj_DifferencesFirstAndSecondVertexFit_deltaTime;   //!
   TBranch        *b_recolns_sigma2_theta;   //!
   TBranch        *b_recolns_phi;   //!
   TBranch        *b_trigger_mask;   //!
   TBranch        *b_gandalf_energy_corrected;   //!
   TBranch        *b_dir_x;   //!
   TBranch        *b_recolns_pos_y;   //!
   TBranch        *b_cosangle_truth_recolns;   //!
   TBranch        *b_recolns_pos_z;   //!
   TBranch        *b_gandalf_spread_pos_y_median;   //!
   TBranch        *b_gandalf_spread_dir_x_iqr;   //!
   TBranch        *b_dusj_llhSinglePMT_sum;   //!
   TBranch        *b__crossValidationTrainingMask_muon_score;   //!
   TBranch        *b_gandalf_spread_dir_z_iqr;   //!
   TBranch        *b_dusj_is_selected_without_anti_noise_cut;   //!
   TBranch        *b_utc_nanoseconds;   //!
   TBranch        *b_gandalf_spread_beta0_median;   //!
   TBranch        *b_gandalf_pos_r;   //!
   TBranch        *b_recolns_pos_xyz;   //!
   TBranch        *b_dusj_llh_sum;   //!
   TBranch        *b_bjorkeny;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_dusj_llhSinglePMT_sum_forNoSignal;   //!
   TBranch        *b_gandalf_spread_chi2_std;   //!
   TBranch        *b_interaction_channel;   //!
   TBranch        *b_gandalf_pos_z;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge;   //!
   TBranch        *b_gandalf_pos_x;   //!
   TBranch        *b_dusj_targetEnergy;   //!
   TBranch        *b_dir_y;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean;   //!
   TBranch        *b_dusj_energyErrorDown_bestLLH;   //!
   TBranch        *b_dir_z;   //!
   TBranch        *b_dusj_EventID;   //!
   TBranch        *b_gandalf_spread_pos_z_mean;   //!
   TBranch        *b_dir_r;   //!
   TBranch        *b_gandalf_spread_n_hits_std;   //!
   TBranch        *b_dusj_best_SecondDusjOrcaVertexFit_FitResult_time;   //!
   TBranch        *b_dusj_shifted30m_pos_y;   //!
   TBranch        *b_dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber;   //!
   TBranch        *b_dusj_shifted30m_pos_z;   //!
   TBranch        *b_pos_xyz;   //!
   TBranch        *b_dusj_zenith;   //!
   TBranch        *b_gandalf_is_good;   //!
   TBranch        *b_dusj_Npmt_maxDeltaT10ns;   //!
   TBranch        *b_gandalf_shifted30m_pos_x;   //!
   TBranch        *b_dusj_dir_xyz;   //!
   TBranch        *b_gandalf_spread_chi2_median;   //!
   TBranch        *b_gandalf_dir_y;   //!
   TBranch        *b_gandalf_dir_x;   //!
   TBranch        *b_gandalf_log10_energy_corrected;   //!
   TBranch        *b_gandalf_dir_z;   //!
   TBranch        *b_distancePointTrack_dusj_truth;   //!
   TBranch        *b_gandalf_loose_is_selected;   //!
   TBranch        *b_gandalf_dir_r;   //!
   TBranch        *b_gandalf_spread_dir_z_median;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean;   //!
   TBranch        *b_length;   //!
   TBranch        *b_gandalf_spread_dir_z_mean;   //!
   TBranch        *b_gandalf_spread_beta1_mean;   //!
   TBranch        *b_gandalf_dir_xyz;   //!
   TBranch        *b_dusj_pos_xyz;   //!
   TBranch        *b_gandalf_spread_pos_y_mean;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct;   //!
   TBranch        *b_recolns_is_selected_without_containment;   //!
   TBranch        *b_dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N;   //!
   TBranch        *b_gandalf_spread_pos_x_mean;   //!
   TBranch        *b_gandalf_cos_theta;   //!
   TBranch        *b_distance_truth_recolns;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_time;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient;   //!
   TBranch        *b_gandalf_spread_pos_z_mad;   //!
   TBranch        *b_gandalf_zenith;   //!
   TBranch        *b_run_id;   //!
   TBranch        *b_dusj_energy_corrected;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY;   //!
   TBranch        *b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_theta;   //!
   TBranch        *b_tracklike_from_Emuon;   //!
   TBranch        *b_dusj_is_selected;   //!
   TBranch        *b_gandalf_spread_chi2_mean;   //!
   TBranch        *b_dusj_DifferencesFirstAndSecondVertexFit_distance;   //!
   TBranch        *b_dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;   //!
   TBranch        *b_dusj_best_FirstDusjOrcaVertexFit_FitResult_time;   //!
   TBranch        *b_dusj_azimuth;   //!
   TBranch        *b_recolns_theta;   //!
   TBranch        *b_recolns_zenith;   //!
   TBranch        *b_gandalf_time;   //!
   TBranch        *b_gandalf_spread_dir_x_std;   //!
   TBranch        *b_gandalf_containment_shifted;   //!
   TBranch        *b_recolns_is_selected;   //!
   TBranch        *b_gandalf_spread_beta1_iqr;   //!
   TBranch        *b_dusj_cosangle_dir_pos;   //!
   TBranch        *b_gandalf_spread_dir_x_median;   //!
   TBranch        *b_dusj_MuonSuppression_decision;   //!
   TBranch        *b_recolns_sigma2_phi;   //!
   TBranch        *b_gandalf_jstart_npe_mip_total;   //!
   TBranch        *b_type;   //!
   TBranch        *b_dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;   //!
   TBranch        *b_dusj_L0AroundL1HitSelection_weight;   //!
   TBranch        *b_recolns_dir_x;   //!
   TBranch        *b_recolns_dir_y;   //!
   TBranch        *b_log_distance4Dshortest_dusj_gandalf;   //!
   TBranch        *b_recolns_dir_r;   //!
   TBranch        *b_cosangle_dusj_recolns;   //!
   TBranch        *b_dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber;   //!
   TBranch        *b_gandalf_spread_beta0_std;   //!
   TBranch        *b_gandalf_spread_n_hits_median;   //!
   TBranch        *b_trigger_counter;   //!
   TBranch        *b_recolns_log10_energy_neutrino;   //!
   TBranch        *b_distancePointTrack_gandalf_truth;   //!
   TBranch        *b_dusj_time;   //!
   TBranch        *b_gandalf_spread_pos_y_mad;   //!
   TBranch        *b_gandalf_chi2;   //!
   TBranch        *b_cos_zenith;   //!
   TBranch        *b_gandalf_jstart_npe_mip;   //!
   TBranch        *b_gandalf_azimuth;   //!
   TBranch        *b_gandalf_spread_dir_z_std;   //!
   TBranch        *b_recolns_energy_neutrino;   //!
   TBranch        *b_dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_z;   //!
   TBranch        *b_track_score;   //!
   TBranch        *b_dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_x;   //!
   TBranch        *b_dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_y;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_z;   //!
   TBranch        *b_mc_t;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_x;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_y;   //!
   TBranch        *b_muon_score;   //!
   TBranch        *b_dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;   //!
   TBranch        *b_gandalf_theta;   //!
   TBranch        *b_dusj_dir_r;   //!
   TBranch        *b_gandalf_n_iter;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean;   //!
   TBranch        *b_gandalf_spread_n_hits_iqr;   //!
   TBranch        *b_recolns_energy_muon;   //!
   TBranch        *b_recolns_containment;   //!
   TBranch        *b_weight_one_year;   //!
   TBranch        *b_distance_dusj_gandalf;   //!
   TBranch        *b_det_id;   //!
   TBranch        *b_recolns_dir_z;   //!
   TBranch        *b_cosangle_truth_dusj;   //!
   TBranch        *b_distance_gandalf_recolns;   //!
   TBranch        *b_zenith;   //!
   TBranch        *b_gandalf_pos_xyz;   //!
   TBranch        *b_dusj_sumMeasuredOMhits;   //!
   TBranch        *b_utc_seconds;   //!
   TBranch        *b_gandalf_spread_pos_y_std;   //!
   TBranch        *b_gandalf_spread_pos_x_std;   //!
   TBranch        *b_dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N;   //!
   TBranch        *b_dusj_sumExpOMhits;   //!
   TBranch        *b_Erange_min;   //!
   TBranch        *b_dusj_pos_r;   //!
   TBranch        *b_gandalf_cos_zenith;   //!
   TBranch        *b_gandalf_spread_beta1_median;   //!
   TBranch        *b_overlays;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits;   //!
   TBranch        *b_id;   //!
   TBranch        *b_dusj_pos_x;   //!
   TBranch        *b_distance_truth_dusj;   //!
   TBranch        *b_dusj_energyErrorUp_bestLLH;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber;   //!
   TBranch        *b_event_id;   //!
   TBranch        *b_dusj_pos_z;   //!
   TBranch        *b_cosangle_truth_gandalf;   //!
   TBranch        *b_dusj_is_good;   //!
   TBranch        *b_gandalf_spread_pos_x_mad;   //!
   TBranch        *b_Erange_max;   //!
   TBranch        *b_gandalf_spread_pos_x_iqr;   //!
   TBranch        *b_gandalf_spread_dir_y_iqr;   //!
   TBranch        *b_gandalf_spread_beta1_std;   //!
   TBranch        *b_dusj_cos_zenith;   //!
   TBranch        *b_gandalf_spread_beta0_mean;   //!
   TBranch        *b_recolns_pos_x;   //!
   TBranch        *b_dusj_pos_y;   //!
   TBranch        *b_dusj_energy_bestLLH;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_dusj_llhBestSinglePMTperDOM_sum;   //!
   TBranch        *b_weight_w3;   //!
   TBranch        *b_weight_w2;   //!
   TBranch        *b_weight_w1;   //!
   TBranch        *b_dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct;   //!
   TBranch        *b_azimuth;   //!
   TBranch        *b_recolns_sin_theta;   //!
   TBranch        *b_gandalf_spread_beta1_mad;   //!
   TBranch        *b_gandalf_spread_chi2_mad;   //!
   TBranch        *b_gandalf_jenergy_chi2;   //!
   TBranch        *b_gandalf_spread_chi2_iqr;   //!
   TBranch        *b_recolns_containment_shifted;   //!
   TBranch        *b_recolns_bjorken_y;   //!
   TBranch        *b_recolns_cos_theta;   //!
   TBranch        *b_recolns_beta;   //!
   TBranch        *b_gandalf_spread_pos_z_iqr;   //!
   TBranch        *b_dusj_cos_theta;   //!
   TBranch        *b_dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber;   //!
   TBranch        *b_gandalf_shifted30m_pos_y;   //!
   TBranch        *b_frame_index;   //!
   TBranch        *b_gandalf_shifted30m_pos_z;   //!
   TBranch        *b_dusj_sumExpFromPoissonOMhits;   //!
   TBranch        *b_gandalf_spread_dir_x_mad;   //!
   TBranch        *b_dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75;   //!
   TBranch        *b__crossValidationTrainingMask_track_score;   //!
   TBranch        *b_gandalf_spread_dir_y_mean;   //!
   TBranch        *b_gandalf_spread_beta0_mad;   //!
   TBranch        *b_dusj_multiplicity;   //!
   TBranch        *b_gandalf_spread_n_hits_mad;   //!
   TBranch        *b_log_distance_dusj_gandalfForDusjTime;   //!
   TBranch        *b_dusj_MuonSuppression_deltaTresQ20Q80;   //!
   TBranch        *b_gandalf_spread_dir_y_median;   //!
   TBranch        *b_time;   //!
   TBranch        *b_recolns_azimuth;   //!
   TBranch        *b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference;   //!
   TBranch        *b_gandalf_spread_dir_y_mad;   //!
   TBranch        *b_recolns_cos_zenith;   //!
   TBranch        *b_dusj_dir_y;   //!
   TBranch        *b_dusj_dir_x;   //!
   TBranch        *b_dusj_dir_z;   //!

   ECAP180401(TTree *tree=0);
   virtual ~ECAP180401();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ECAP180401_cxx
ECAP180401::ECAP180401(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("irodsdata/pid_output_atm_neutrino_atm_muon_pure_noise_shiftedVertexEventSelection_180401.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("irodsdata/pid_output_atm_neutrino_atm_muon_pure_noise_shiftedVertexEventSelection_180401.root");
      }
      f->GetObject("PID",tree);

   }
   Init(tree);
}

ECAP180401::~ECAP180401()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ECAP180401::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ECAP180401::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ECAP180401::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("distance_dusj_recolns", &distance_dusj_recolns, &b_distance_dusj_recolns);
   fChain->SetBranchAddress("recolns_pos_r", &recolns_pos_r, &b_recolns_pos_r);
   fChain->SetBranchAddress("dusj_theta", &dusj_theta, &b_dusj_theta);
   fChain->SetBranchAddress("recolns_dir_xyz", &recolns_dir_xyz, &b_recolns_dir_xyz);
   fChain->SetBranchAddress("recolns_shifted30m_pos_z", &recolns_shifted30m_pos_z, &b_recolns_shifted30m_pos_z);
   fChain->SetBranchAddress("dusj_log10_energy_corrected", &dusj_log10_energy_corrected, &b_dusj_log10_energy_corrected);
   fChain->SetBranchAddress("recolns_n_hits_used", &recolns_n_hits_used, &b_recolns_n_hits_used);
   fChain->SetBranchAddress("gandalf_phi", &gandalf_phi, &b_gandalf_phi);
   fChain->SetBranchAddress("gandalf_is_selected_without_containment", &gandalf_is_selected_without_containment, &b_gandalf_is_selected_without_containment);
   fChain->SetBranchAddress("dusj_sumMeasuredPMThits", &dusj_sumMeasuredPMThits, &b_dusj_sumMeasuredPMThits);
   fChain->SetBranchAddress("dusj_MuonSuppression_enoughHits", &dusj_MuonSuppression_enoughHits, &b_dusj_MuonSuppression_enoughHits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_energy", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_energy, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_energy);
   fChain->SetBranchAddress("recolns_shifted30m_pos_y", &recolns_shifted30m_pos_y, &b_recolns_shifted30m_pos_y);
   fChain->SetBranchAddress("recolns_quality", &recolns_quality, &b_recolns_quality);
   fChain->SetBranchAddress("_crossValidationTrainingMask_noise_score", &_crossValidationTrainingMask_noise_score, &b__crossValidationTrainingMask_noise_score);
   fChain->SetBranchAddress("distance_truth_gandalf", &distance_truth_gandalf, &b_distance_truth_gandalf);
   fChain->SetBranchAddress("dusj_L0AroundL1HitSelection__weight_withoutSelfSquared", &dusj_L0AroundL1HitSelection__weight_withoutSelfSquared, &b_dusj_L0AroundL1HitSelection__weight_withoutSelfSquared);
   fChain->SetBranchAddress("dusj_shifted30m_pos_x", &dusj_shifted30m_pos_x, &b_dusj_shifted30m_pos_x);
   fChain->SetBranchAddress("dusj_passes_anti_noise_cut", &dusj_passes_anti_noise_cut, &b_dusj_passes_anti_noise_cut);
   fChain->SetBranchAddress("dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80", &dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80, &b_dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80);
   fChain->SetBranchAddress("dusj_phi", &dusj_phi, &b_dusj_phi);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge);
   fChain->SetBranchAddress("gandalf_spread_pos_x_median", &gandalf_spread_pos_x_median, &b_gandalf_spread_pos_x_median);
   fChain->SetBranchAddress("mc_id", &mc_id, &b_mc_id);
   fChain->SetBranchAddress("dusj_L0AroundL1HitSelection_weight_withoutSelfSquared", &dusj_L0AroundL1HitSelection_weight_withoutSelfSquared, &b_dusj_L0AroundL1HitSelection_weight_withoutSelfSquared);
   fChain->SetBranchAddress("gandalf_spread_pos_y_iqr", &gandalf_spread_pos_y_iqr, &b_gandalf_spread_pos_y_iqr);
   fChain->SetBranchAddress("dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference", &dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference, &b_dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference);
   fChain->SetBranchAddress("dusj_llhBestSinglePMTperDOM__sum_forNoSignal", &dusj_llhBestSinglePMTperDOM__sum_forNoSignal, &b_dusj_llhBestSinglePMTperDOM__sum_forNoSignal);
   fChain->SetBranchAddress("dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_y", &dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_y, &b_dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_y);
   fChain->SetBranchAddress("dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_x", &dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_x, &b_dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_x);
   fChain->SetBranchAddress("dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_z", &dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_z, &b_dusj_best_FirstDusjOrcaVertexFit_FitResult_pos_z);
   fChain->SetBranchAddress("log10_energy", &log10_energy, &b_log10_energy);
   fChain->SetBranchAddress("gandalf_spread_dir_z_mad", &gandalf_spread_dir_z_mad, &b_gandalf_spread_dir_z_mad);
   fChain->SetBranchAddress("gandalf_is_selected", &gandalf_is_selected, &b_gandalf_is_selected);
   fChain->SetBranchAddress("dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N", &dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N, &b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N);
   fChain->SetBranchAddress("noise_score", &noise_score, &b_noise_score);
   fChain->SetBranchAddress("gandalf_jstart_length", &gandalf_jstart_length, &b_gandalf_jstart_length);
   fChain->SetBranchAddress("dusj_Nom", &dusj_Nom, &b_dusj_Nom);
   fChain->SetBranchAddress("cosangle_gandalf_recolns", &cosangle_gandalf_recolns, &b_cosangle_gandalf_recolns);
   fChain->SetBranchAddress("gandalf_cosangle_dir_pos", &gandalf_cosangle_dir_pos, &b_gandalf_cosangle_dir_pos);
   fChain->SetBranchAddress("gandalf_lambda", &gandalf_lambda, &b_gandalf_lambda);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_phi", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_phi, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_phi);
   fChain->SetBranchAddress("selected_in_ECAP_PID", &selected_in_ECAP_PID, &b_selected_in_ECAP_PID);
   fChain->SetBranchAddress("dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber", &dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber, &b_dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber);
   fChain->SetBranchAddress("dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference", &dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference, &b_dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference);
   fChain->SetBranchAddress("gandalf_containment", &gandalf_containment, &b_gandalf_containment);
   fChain->SetBranchAddress("gandalf_spread_dir_x_mean", &gandalf_spread_dir_x_mean, &b_gandalf_spread_dir_x_mean);
   fChain->SetBranchAddress("pos_z", &pos_z, &b_pos_z);
   fChain->SetBranchAddress("pos_x", &pos_x, &b_pos_x);
   fChain->SetBranchAddress("pos_y", &pos_y, &b_pos_y);
   fChain->SetBranchAddress("gandalf_beta1", &gandalf_beta1, &b_gandalf_beta1);
   fChain->SetBranchAddress("dusj_energy", &dusj_energy, &b_dusj_energy);
   fChain->SetBranchAddress("gandalf_beta0", &gandalf_beta0, &b_gandalf_beta0);
   fChain->SetBranchAddress("pos_r", &pos_r, &b_pos_r);
   fChain->SetBranchAddress("dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75", &dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75, &b_dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75);
   fChain->SetBranchAddress("gandalf_spread_dir_y_std", &gandalf_spread_dir_y_std, &b_gandalf_spread_dir_y_std);
   fChain->SetBranchAddress("recolns_loose_is_selected", &recolns_loose_is_selected, &b_recolns_loose_is_selected);
   fChain->SetBranchAddress("recolns_cosangle_dir_pos", &recolns_cosangle_dir_pos, &b_recolns_cosangle_dir_pos);
   fChain->SetBranchAddress("dusj_Fork_muonSuppression_decision", &dusj_Fork_muonSuppression_decision, &b_dusj_Fork_muonSuppression_decision);
   fChain->SetBranchAddress("gandalf_spread_pos_z_std", &gandalf_spread_pos_z_std, &b_gandalf_spread_pos_z_std);
   fChain->SetBranchAddress("gandalf_spread_pos_z_median", &gandalf_spread_pos_z_median, &b_gandalf_spread_pos_z_median);
   fChain->SetBranchAddress("cosangle_dusj_gandalf", &cosangle_dusj_gandalf, &b_cosangle_dusj_gandalf);
   fChain->SetBranchAddress("gandalf_jenergy_energy", &gandalf_jenergy_energy, &b_gandalf_jenergy_energy);
   fChain->SetBranchAddress("livetime_sec", &livetime_sec, &b_livetime_sec);
   fChain->SetBranchAddress("recolns_is_good", &recolns_is_good, &b_recolns_is_good);
   fChain->SetBranchAddress("cos_theta", &cos_theta, &b_cos_theta);
   fChain->SetBranchAddress("is_cc", &is_cc, &b_is_cc);
   fChain->SetBranchAddress("gandalf_n_hits", &gandalf_n_hits, &b_gandalf_n_hits);
   fChain->SetBranchAddress("n_files_gen", &n_files_gen, &b_n_files_gen);
   fChain->SetBranchAddress("gandalf_spread_beta0_iqr", &gandalf_spread_beta0_iqr, &b_gandalf_spread_beta0_iqr);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("dusj_Npmt", &dusj_Npmt, &b_dusj_Npmt);
   fChain->SetBranchAddress("Emuon_from_neutrino_interaction", &Emuon_from_neutrino_interaction, &b_Emuon_from_neutrino_interaction);
   fChain->SetBranchAddress("recolns_shifted30m_pos_x", &recolns_shifted30m_pos_x, &b_recolns_shifted30m_pos_x);
   fChain->SetBranchAddress("recolns_log10_energy_muon", &recolns_log10_energy_muon, &b_recolns_log10_energy_muon);
   fChain->SetBranchAddress("is_neutrino", &is_neutrino, &b_is_neutrino);
   fChain->SetBranchAddress("gandalf_pos_y", &gandalf_pos_y, &b_gandalf_pos_y);
   fChain->SetBranchAddress("dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult", &dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult, &b_dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult);
   fChain->SetBranchAddress("n_events_gen", &n_events_gen, &b_n_events_gen);
   fChain->SetBranchAddress("dir_xyz", &dir_xyz, &b_dir_xyz);
   fChain->SetBranchAddress("gandalf_spread_n_hits_mean", &gandalf_spread_n_hits_mean, &b_gandalf_spread_n_hits_mean);
   fChain->SetBranchAddress("dusj_DifferencesFirstAndSecondVertexFit_deltaTime", &dusj_DifferencesFirstAndSecondVertexFit_deltaTime, &b_dusj_DifferencesFirstAndSecondVertexFit_deltaTime);
   fChain->SetBranchAddress("recolns_sigma2_theta", &recolns_sigma2_theta, &b_recolns_sigma2_theta);
   fChain->SetBranchAddress("recolns_phi", &recolns_phi, &b_recolns_phi);
   fChain->SetBranchAddress("trigger_mask", &trigger_mask, &b_trigger_mask);
   fChain->SetBranchAddress("gandalf_energy_corrected", &gandalf_energy_corrected, &b_gandalf_energy_corrected);
   fChain->SetBranchAddress("dir_x", &dir_x, &b_dir_x);
   fChain->SetBranchAddress("recolns_pos_y", &recolns_pos_y, &b_recolns_pos_y);
   fChain->SetBranchAddress("cosangle_truth_recolns", &cosangle_truth_recolns, &b_cosangle_truth_recolns);
   fChain->SetBranchAddress("recolns_pos_z", &recolns_pos_z, &b_recolns_pos_z);
   fChain->SetBranchAddress("gandalf_spread_pos_y_median", &gandalf_spread_pos_y_median, &b_gandalf_spread_pos_y_median);
   fChain->SetBranchAddress("gandalf_spread_dir_x_iqr", &gandalf_spread_dir_x_iqr, &b_gandalf_spread_dir_x_iqr);
   fChain->SetBranchAddress("dusj_llhSinglePMT_sum", &dusj_llhSinglePMT_sum, &b_dusj_llhSinglePMT_sum);
   fChain->SetBranchAddress("_crossValidationTrainingMask_muon_score", &_crossValidationTrainingMask_muon_score, &b__crossValidationTrainingMask_muon_score);
   fChain->SetBranchAddress("gandalf_spread_dir_z_iqr", &gandalf_spread_dir_z_iqr, &b_gandalf_spread_dir_z_iqr);
   fChain->SetBranchAddress("dusj_is_selected_without_anti_noise_cut", &dusj_is_selected_without_anti_noise_cut, &b_dusj_is_selected_without_anti_noise_cut);
   fChain->SetBranchAddress("utc_nanoseconds", &utc_nanoseconds, &b_utc_nanoseconds);
   fChain->SetBranchAddress("gandalf_spread_beta0_median", &gandalf_spread_beta0_median, &b_gandalf_spread_beta0_median);
   fChain->SetBranchAddress("gandalf_pos_r", &gandalf_pos_r, &b_gandalf_pos_r);
   fChain->SetBranchAddress("recolns_pos_xyz", &recolns_pos_xyz, &b_recolns_pos_xyz);
   fChain->SetBranchAddress("dusj_llh_sum", &dusj_llh_sum, &b_dusj_llh_sum);
   fChain->SetBranchAddress("bjorkeny", &bjorkeny, &b_bjorkeny);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("dusj_llhSinglePMT_sum_forNoSignal", &dusj_llhSinglePMT_sum_forNoSignal, &b_dusj_llhSinglePMT_sum_forNoSignal);
   fChain->SetBranchAddress("gandalf_spread_chi2_std", &gandalf_spread_chi2_std, &b_gandalf_spread_chi2_std);
   fChain->SetBranchAddress("interaction_channel", &interaction_channel, &b_interaction_channel);
   fChain->SetBranchAddress("gandalf_pos_z", &gandalf_pos_z, &b_gandalf_pos_z);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge);
   fChain->SetBranchAddress("gandalf_pos_x", &gandalf_pos_x, &b_gandalf_pos_x);
   fChain->SetBranchAddress("dusj_targetEnergy", &dusj_targetEnergy, &b_dusj_targetEnergy);
   fChain->SetBranchAddress("dir_y", &dir_y, &b_dir_y);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean);
   fChain->SetBranchAddress("dusj_energyErrorDown_bestLLH", &dusj_energyErrorDown_bestLLH, &b_dusj_energyErrorDown_bestLLH);
   fChain->SetBranchAddress("dir_z", &dir_z, &b_dir_z);
   fChain->SetBranchAddress("dusj_EventID", &dusj_EventID, &b_dusj_EventID);
   fChain->SetBranchAddress("gandalf_spread_pos_z_mean", &gandalf_spread_pos_z_mean, &b_gandalf_spread_pos_z_mean);
   fChain->SetBranchAddress("dir_r", &dir_r, &b_dir_r);
   fChain->SetBranchAddress("gandalf_spread_n_hits_std", &gandalf_spread_n_hits_std, &b_gandalf_spread_n_hits_std);
   fChain->SetBranchAddress("dusj_best_SecondDusjOrcaVertexFit_FitResult_time", &dusj_best_SecondDusjOrcaVertexFit_FitResult_time, &b_dusj_best_SecondDusjOrcaVertexFit_FitResult_time);
   fChain->SetBranchAddress("dusj_shifted30m_pos_y", &dusj_shifted30m_pos_y, &b_dusj_shifted30m_pos_y);
   fChain->SetBranchAddress("dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber", &dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber, &b_dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber);
   fChain->SetBranchAddress("dusj_shifted30m_pos_z", &dusj_shifted30m_pos_z, &b_dusj_shifted30m_pos_z);
   fChain->SetBranchAddress("pos_xyz", &pos_xyz, &b_pos_xyz);
   fChain->SetBranchAddress("dusj_zenith", &dusj_zenith, &b_dusj_zenith);
   fChain->SetBranchAddress("gandalf_is_good", &gandalf_is_good, &b_gandalf_is_good);
   fChain->SetBranchAddress("dusj_Npmt_maxDeltaT10ns", &dusj_Npmt_maxDeltaT10ns, &b_dusj_Npmt_maxDeltaT10ns);
   fChain->SetBranchAddress("gandalf_shifted30m_pos_x", &gandalf_shifted30m_pos_x, &b_gandalf_shifted30m_pos_x);
   fChain->SetBranchAddress("dusj_dir_xyz", &dusj_dir_xyz, &b_dusj_dir_xyz);
   fChain->SetBranchAddress("gandalf_spread_chi2_median", &gandalf_spread_chi2_median, &b_gandalf_spread_chi2_median);
   fChain->SetBranchAddress("gandalf_dir_y", &gandalf_dir_y, &b_gandalf_dir_y);
   fChain->SetBranchAddress("gandalf_dir_x", &gandalf_dir_x, &b_gandalf_dir_x);
   fChain->SetBranchAddress("gandalf_log10_energy_corrected", &gandalf_log10_energy_corrected, &b_gandalf_log10_energy_corrected);
   fChain->SetBranchAddress("gandalf_dir_z", &gandalf_dir_z, &b_gandalf_dir_z);
   fChain->SetBranchAddress("distancePointTrack_dusj_truth", &distancePointTrack_dusj_truth, &b_distancePointTrack_dusj_truth);
   fChain->SetBranchAddress("gandalf_loose_is_selected", &gandalf_loose_is_selected, &b_gandalf_loose_is_selected);
   fChain->SetBranchAddress("gandalf_dir_r", &gandalf_dir_r, &b_gandalf_dir_r);
   fChain->SetBranchAddress("gandalf_spread_dir_z_median", &gandalf_spread_dir_z_median, &b_gandalf_spread_dir_z_median);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean);
   fChain->SetBranchAddress("length", &length, &b_length);
   fChain->SetBranchAddress("gandalf_spread_dir_z_mean", &gandalf_spread_dir_z_mean, &b_gandalf_spread_dir_z_mean);
   fChain->SetBranchAddress("gandalf_spread_beta1_mean", &gandalf_spread_beta1_mean, &b_gandalf_spread_beta1_mean);
   fChain->SetBranchAddress("gandalf_dir_xyz", &gandalf_dir_xyz, &b_gandalf_dir_xyz);
   fChain->SetBranchAddress("dusj_pos_xyz", &dusj_pos_xyz, &b_dusj_pos_xyz);
   fChain->SetBranchAddress("gandalf_spread_pos_y_mean", &gandalf_spread_pos_y_mean, &b_gandalf_spread_pos_y_mean);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct);
   fChain->SetBranchAddress("recolns_is_selected_without_containment", &recolns_is_selected_without_containment, &b_recolns_is_selected_without_containment);
   fChain->SetBranchAddress("dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N", &dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N, &b_dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N);
   fChain->SetBranchAddress("gandalf_spread_pos_x_mean", &gandalf_spread_pos_x_mean, &b_gandalf_spread_pos_x_mean);
   fChain->SetBranchAddress("gandalf_cos_theta", &gandalf_cos_theta, &b_gandalf_cos_theta);
   fChain->SetBranchAddress("distance_truth_recolns", &distance_truth_recolns, &b_distance_truth_recolns);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_time", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_time, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_time);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient);
   fChain->SetBranchAddress("gandalf_spread_pos_z_mad", &gandalf_spread_pos_z_mad, &b_gandalf_spread_pos_z_mad);
   fChain->SetBranchAddress("gandalf_zenith", &gandalf_zenith, &b_gandalf_zenith);
   fChain->SetBranchAddress("run_id", &run_id, &b_run_id);
   fChain->SetBranchAddress("dusj_energy_corrected", &dusj_energy_corrected, &b_dusj_energy_corrected);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY);
   fChain->SetBranchAddress("dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference", &dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference, &b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_theta", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_theta, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_theta);
   fChain->SetBranchAddress("tracklike_from_Emuon", &tracklike_from_Emuon, &b_tracklike_from_Emuon);
   fChain->SetBranchAddress("dusj_is_selected", &dusj_is_selected, &b_dusj_is_selected);
   fChain->SetBranchAddress("gandalf_spread_chi2_mean", &gandalf_spread_chi2_mean, &b_gandalf_spread_chi2_mean);
   fChain->SetBranchAddress("dusj_DifferencesFirstAndSecondVertexFit_distance", &dusj_DifferencesFirstAndSecondVertexFit_distance, &b_dusj_DifferencesFirstAndSecondVertexFit_distance);
   fChain->SetBranchAddress("dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult", &dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult, &b_dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult);
   fChain->SetBranchAddress("dusj_best_FirstDusjOrcaVertexFit_FitResult_time", &dusj_best_FirstDusjOrcaVertexFit_FitResult_time, &b_dusj_best_FirstDusjOrcaVertexFit_FitResult_time);
   fChain->SetBranchAddress("dusj_azimuth", &dusj_azimuth, &b_dusj_azimuth);
   fChain->SetBranchAddress("recolns_theta", &recolns_theta, &b_recolns_theta);
   fChain->SetBranchAddress("recolns_zenith", &recolns_zenith, &b_recolns_zenith);
   fChain->SetBranchAddress("gandalf_time", &gandalf_time, &b_gandalf_time);
   fChain->SetBranchAddress("gandalf_spread_dir_x_std", &gandalf_spread_dir_x_std, &b_gandalf_spread_dir_x_std);
   fChain->SetBranchAddress("gandalf_containment_shifted", &gandalf_containment_shifted, &b_gandalf_containment_shifted);
   fChain->SetBranchAddress("recolns_is_selected", &recolns_is_selected, &b_recolns_is_selected);
   fChain->SetBranchAddress("gandalf_spread_beta1_iqr", &gandalf_spread_beta1_iqr, &b_gandalf_spread_beta1_iqr);
   fChain->SetBranchAddress("dusj_cosangle_dir_pos", &dusj_cosangle_dir_pos, &b_dusj_cosangle_dir_pos);
   fChain->SetBranchAddress("gandalf_spread_dir_x_median", &gandalf_spread_dir_x_median, &b_gandalf_spread_dir_x_median);
   fChain->SetBranchAddress("dusj_MuonSuppression_decision", &dusj_MuonSuppression_decision, &b_dusj_MuonSuppression_decision);
   fChain->SetBranchAddress("recolns_sigma2_phi", &recolns_sigma2_phi, &b_recolns_sigma2_phi);
   fChain->SetBranchAddress("gandalf_jstart_npe_mip_total", &gandalf_jstart_npe_mip_total, &b_gandalf_jstart_npe_mip_total);
   fChain->SetBranchAddress("type", &type, &b_type);
   fChain->SetBranchAddress("dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult", &dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult, &b_dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult);
   fChain->SetBranchAddress("dusj_L0AroundL1HitSelection_weight", &dusj_L0AroundL1HitSelection_weight, &b_dusj_L0AroundL1HitSelection_weight);
   fChain->SetBranchAddress("recolns_dir_x", &recolns_dir_x, &b_recolns_dir_x);
   fChain->SetBranchAddress("recolns_dir_y", &recolns_dir_y, &b_recolns_dir_y);
   fChain->SetBranchAddress("log_distance4Dshortest_dusj_gandalf", &log_distance4Dshortest_dusj_gandalf, &b_log_distance4Dshortest_dusj_gandalf);
   fChain->SetBranchAddress("recolns_dir_r", &recolns_dir_r, &b_recolns_dir_r);
   fChain->SetBranchAddress("cosangle_dusj_recolns", &cosangle_dusj_recolns, &b_cosangle_dusj_recolns);
   fChain->SetBranchAddress("dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber", &dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber, &b_dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber);
   fChain->SetBranchAddress("gandalf_spread_beta0_std", &gandalf_spread_beta0_std, &b_gandalf_spread_beta0_std);
   fChain->SetBranchAddress("gandalf_spread_n_hits_median", &gandalf_spread_n_hits_median, &b_gandalf_spread_n_hits_median);
   fChain->SetBranchAddress("trigger_counter", &trigger_counter, &b_trigger_counter);
   fChain->SetBranchAddress("recolns_log10_energy_neutrino", &recolns_log10_energy_neutrino, &b_recolns_log10_energy_neutrino);
   fChain->SetBranchAddress("distancePointTrack_gandalf_truth", &distancePointTrack_gandalf_truth, &b_distancePointTrack_gandalf_truth);
   fChain->SetBranchAddress("dusj_time", &dusj_time, &b_dusj_time);
   fChain->SetBranchAddress("gandalf_spread_pos_y_mad", &gandalf_spread_pos_y_mad, &b_gandalf_spread_pos_y_mad);
   fChain->SetBranchAddress("gandalf_chi2", &gandalf_chi2, &b_gandalf_chi2);
   fChain->SetBranchAddress("cos_zenith", &cos_zenith, &b_cos_zenith);
   fChain->SetBranchAddress("gandalf_jstart_npe_mip", &gandalf_jstart_npe_mip, &b_gandalf_jstart_npe_mip);
   fChain->SetBranchAddress("gandalf_azimuth", &gandalf_azimuth, &b_gandalf_azimuth);
   fChain->SetBranchAddress("gandalf_spread_dir_z_std", &gandalf_spread_dir_z_std, &b_gandalf_spread_dir_z_std);
   fChain->SetBranchAddress("recolns_energy_neutrino", &recolns_energy_neutrino, &b_recolns_energy_neutrino);
   fChain->SetBranchAddress("dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_z", &dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_z, &b_dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_z);
   fChain->SetBranchAddress("track_score", &track_score, &b_track_score);
   fChain->SetBranchAddress("dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_x", &dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_x, &b_dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_x);
   fChain->SetBranchAddress("dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_y", &dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_y, &b_dusj_best_SecondDusjOrcaVertexFit_FitResult_pos_y);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_z", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_z, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_z);
   fChain->SetBranchAddress("mc_t", &mc_t, &b_mc_t);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_x", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_x, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_x);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_y", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_y, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_pos_y);
   fChain->SetBranchAddress("muon_score", &muon_score, &b_muon_score);
   fChain->SetBranchAddress("dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult", &dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult, &b_dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult);
   fChain->SetBranchAddress("gandalf_theta", &gandalf_theta, &b_gandalf_theta);
   fChain->SetBranchAddress("dusj_dir_r", &dusj_dir_r, &b_dusj_dir_r);
   fChain->SetBranchAddress("gandalf_n_iter", &gandalf_n_iter, &b_gandalf_n_iter);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean);
   fChain->SetBranchAddress("gandalf_spread_n_hits_iqr", &gandalf_spread_n_hits_iqr, &b_gandalf_spread_n_hits_iqr);
   fChain->SetBranchAddress("recolns_energy_muon", &recolns_energy_muon, &b_recolns_energy_muon);
   fChain->SetBranchAddress("recolns_containment", &recolns_containment, &b_recolns_containment);
   fChain->SetBranchAddress("weight_one_year", &weight_one_year, &b_weight_one_year);
   fChain->SetBranchAddress("distance_dusj_gandalf", &distance_dusj_gandalf, &b_distance_dusj_gandalf);
   fChain->SetBranchAddress("det_id", &det_id, &b_det_id);
   fChain->SetBranchAddress("recolns_dir_z", &recolns_dir_z, &b_recolns_dir_z);
   fChain->SetBranchAddress("cosangle_truth_dusj", &cosangle_truth_dusj, &b_cosangle_truth_dusj);
   fChain->SetBranchAddress("distance_gandalf_recolns", &distance_gandalf_recolns, &b_distance_gandalf_recolns);
   fChain->SetBranchAddress("zenith", &zenith, &b_zenith);
   fChain->SetBranchAddress("gandalf_pos_xyz", &gandalf_pos_xyz, &b_gandalf_pos_xyz);
   fChain->SetBranchAddress("dusj_sumMeasuredOMhits", &dusj_sumMeasuredOMhits, &b_dusj_sumMeasuredOMhits);
   fChain->SetBranchAddress("utc_seconds", &utc_seconds, &b_utc_seconds);
   fChain->SetBranchAddress("gandalf_spread_pos_y_std", &gandalf_spread_pos_y_std, &b_gandalf_spread_pos_y_std);
   fChain->SetBranchAddress("gandalf_spread_pos_x_std", &gandalf_spread_pos_x_std, &b_gandalf_spread_pos_x_std);
   fChain->SetBranchAddress("dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N", &dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N, &b_dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N);
   fChain->SetBranchAddress("dusj_sumExpOMhits", &dusj_sumExpOMhits, &b_dusj_sumExpOMhits);
   fChain->SetBranchAddress("Erange_min", &Erange_min, &b_Erange_min);
   fChain->SetBranchAddress("dusj_pos_r", &dusj_pos_r, &b_dusj_pos_r);
   fChain->SetBranchAddress("gandalf_cos_zenith", &gandalf_cos_zenith, &b_gandalf_cos_zenith);
   fChain->SetBranchAddress("gandalf_spread_beta1_median", &gandalf_spread_beta1_median, &b_gandalf_spread_beta1_median);
   fChain->SetBranchAddress("overlays", &overlays, &b_overlays);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits);
   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("dusj_pos_x", &dusj_pos_x, &b_dusj_pos_x);
   fChain->SetBranchAddress("distance_truth_dusj", &distance_truth_dusj, &b_distance_truth_dusj);
   fChain->SetBranchAddress("dusj_energyErrorUp_bestLLH", &dusj_energyErrorUp_bestLLH, &b_dusj_energyErrorUp_bestLLH);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber);
   fChain->SetBranchAddress("event_id", &event_id, &b_event_id);
   fChain->SetBranchAddress("dusj_pos_z", &dusj_pos_z, &b_dusj_pos_z);
   fChain->SetBranchAddress("cosangle_truth_gandalf", &cosangle_truth_gandalf, &b_cosangle_truth_gandalf);
   fChain->SetBranchAddress("dusj_is_good", &dusj_is_good, &b_dusj_is_good);
   fChain->SetBranchAddress("gandalf_spread_pos_x_mad", &gandalf_spread_pos_x_mad, &b_gandalf_spread_pos_x_mad);
   fChain->SetBranchAddress("Erange_max", &Erange_max, &b_Erange_max);
   fChain->SetBranchAddress("gandalf_spread_pos_x_iqr", &gandalf_spread_pos_x_iqr, &b_gandalf_spread_pos_x_iqr);
   fChain->SetBranchAddress("gandalf_spread_dir_y_iqr", &gandalf_spread_dir_y_iqr, &b_gandalf_spread_dir_y_iqr);
   fChain->SetBranchAddress("gandalf_spread_beta1_std", &gandalf_spread_beta1_std, &b_gandalf_spread_beta1_std);
   fChain->SetBranchAddress("dusj_cos_zenith", &dusj_cos_zenith, &b_dusj_cos_zenith);
   fChain->SetBranchAddress("gandalf_spread_beta0_mean", &gandalf_spread_beta0_mean, &b_gandalf_spread_beta0_mean);
   fChain->SetBranchAddress("recolns_pos_x", &recolns_pos_x, &b_recolns_pos_x);
   fChain->SetBranchAddress("dusj_pos_y", &dusj_pos_y, &b_dusj_pos_y);
   fChain->SetBranchAddress("dusj_energy_bestLLH", &dusj_energy_bestLLH, &b_dusj_energy_bestLLH);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("dusj_llhBestSinglePMTperDOM_sum", &dusj_llhBestSinglePMTperDOM_sum, &b_dusj_llhBestSinglePMTperDOM_sum);
   fChain->SetBranchAddress("weight_w3", &weight_w3, &b_weight_w3);
   fChain->SetBranchAddress("weight_w2", &weight_w2, &b_weight_w2);
   fChain->SetBranchAddress("weight_w1", &weight_w1, &b_weight_w1);
   fChain->SetBranchAddress("dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80", &dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80, &b_dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct);
   fChain->SetBranchAddress("azimuth", &azimuth, &b_azimuth);
   fChain->SetBranchAddress("recolns_sin_theta", &recolns_sin_theta, &b_recolns_sin_theta);
   fChain->SetBranchAddress("gandalf_spread_beta1_mad", &gandalf_spread_beta1_mad, &b_gandalf_spread_beta1_mad);
   fChain->SetBranchAddress("gandalf_spread_chi2_mad", &gandalf_spread_chi2_mad, &b_gandalf_spread_chi2_mad);
   fChain->SetBranchAddress("gandalf_jenergy_chi2", &gandalf_jenergy_chi2, &b_gandalf_jenergy_chi2);
   fChain->SetBranchAddress("gandalf_spread_chi2_iqr", &gandalf_spread_chi2_iqr, &b_gandalf_spread_chi2_iqr);
   fChain->SetBranchAddress("recolns_containment_shifted", &recolns_containment_shifted, &b_recolns_containment_shifted);
   fChain->SetBranchAddress("recolns_bjorken_y", &recolns_bjorken_y, &b_recolns_bjorken_y);
   fChain->SetBranchAddress("recolns_cos_theta", &recolns_cos_theta, &b_recolns_cos_theta);
   fChain->SetBranchAddress("recolns_beta", &recolns_beta, &b_recolns_beta);
   fChain->SetBranchAddress("gandalf_spread_pos_z_iqr", &gandalf_spread_pos_z_iqr, &b_gandalf_spread_pos_z_iqr);
   fChain->SetBranchAddress("dusj_cos_theta", &dusj_cos_theta, &b_dusj_cos_theta);
   fChain->SetBranchAddress("dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber", &dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber, &b_dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber);
   fChain->SetBranchAddress("gandalf_shifted30m_pos_y", &gandalf_shifted30m_pos_y, &b_gandalf_shifted30m_pos_y);
   fChain->SetBranchAddress("frame_index", &frame_index, &b_frame_index);
   fChain->SetBranchAddress("gandalf_shifted30m_pos_z", &gandalf_shifted30m_pos_z, &b_gandalf_shifted30m_pos_z);
   fChain->SetBranchAddress("dusj_sumExpFromPoissonOMhits", &dusj_sumExpFromPoissonOMhits, &b_dusj_sumExpFromPoissonOMhits);
   fChain->SetBranchAddress("gandalf_spread_dir_x_mad", &gandalf_spread_dir_x_mad, &b_gandalf_spread_dir_x_mad);
   fChain->SetBranchAddress("dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75", &dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75, &b_dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75);
   fChain->SetBranchAddress("_crossValidationTrainingMask_track_score", &_crossValidationTrainingMask_track_score, &b__crossValidationTrainingMask_track_score);
   fChain->SetBranchAddress("gandalf_spread_dir_y_mean", &gandalf_spread_dir_y_mean, &b_gandalf_spread_dir_y_mean);
   fChain->SetBranchAddress("gandalf_spread_beta0_mad", &gandalf_spread_beta0_mad, &b_gandalf_spread_beta0_mad);
   fChain->SetBranchAddress("dusj_multiplicity", &dusj_multiplicity, &b_dusj_multiplicity);
   fChain->SetBranchAddress("gandalf_spread_n_hits_mad", &gandalf_spread_n_hits_mad, &b_gandalf_spread_n_hits_mad);
   fChain->SetBranchAddress("log_distance_dusj_gandalfForDusjTime", &log_distance_dusj_gandalfForDusjTime, &b_log_distance_dusj_gandalfForDusjTime);
   fChain->SetBranchAddress("dusj_MuonSuppression_deltaTresQ20Q80", &dusj_MuonSuppression_deltaTresQ20Q80, &b_dusj_MuonSuppression_deltaTresQ20Q80);
   fChain->SetBranchAddress("gandalf_spread_dir_y_median", &gandalf_spread_dir_y_median, &b_gandalf_spread_dir_y_median);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("recolns_azimuth", &recolns_azimuth, &b_recolns_azimuth);
   fChain->SetBranchAddress("dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference", &dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference, &b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference);
   fChain->SetBranchAddress("gandalf_spread_dir_y_mad", &gandalf_spread_dir_y_mad, &b_gandalf_spread_dir_y_mad);
   fChain->SetBranchAddress("recolns_cos_zenith", &recolns_cos_zenith, &b_recolns_cos_zenith);
   fChain->SetBranchAddress("dusj_dir_y", &dusj_dir_y, &b_dusj_dir_y);
   fChain->SetBranchAddress("dusj_dir_x", &dusj_dir_x, &b_dusj_dir_x);
   fChain->SetBranchAddress("dusj_dir_z", &dusj_dir_z, &b_dusj_dir_z);
   Notify();
}

Bool_t ECAP180401::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ECAP180401::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ECAP180401::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ECAP180401_cxx
