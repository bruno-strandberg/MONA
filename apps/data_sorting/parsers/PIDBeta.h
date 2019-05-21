//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 30 10:52:34 2019 by ROOT version 6.10/02
// from TTree PID/PID
// found on file: ../../data/pid_result_13Oct2018_ORCA7.root
//////////////////////////////////////////////////////////

#ifndef PIDBeta_h
#define PIDBeta_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class PIDBeta {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        gandalf_JENERGY_CHI2;
   Double_t        gandalf_JENERGY_ENERGY;
   Double_t        gandalf_JENERGY_MUON_RANGE_METRES;
   Double_t        gandalf_JENERGY_NDF;
   Double_t        gandalf_JENERGY_NOISE_LIKELIHOOD;
   Double_t        gandalf_JENERGY_NUMBER_OF_HITS;
   Double_t        gandalf_JGANDALF_BETA0_RAD;
   Double_t        gandalf_JGANDALF_BETA1_RAD;
   Double_t        gandalf_JGANDALF_CHI2;
   Double_t        gandalf_JGANDALF_LAMBDA;
   Double_t        gandalf_JGANDALF_NUMBER_OF_HITS;
   Double_t        gandalf_JGANDALF_NUMBER_OF_ITERATIONS;
   Double_t        gandalf_JSTART_LENGTH_METRES;
   Double_t        gandalf_JSTART_NPE_MIP;
   Double_t        gandalf_JSTART_NPE_MIP_TOTAL;
   Double_t        gandalf_JVETO_NPE;
   Double_t        gandalf_JVETO_NUMBER_OF_HITS;
   Double_t        gandalf_cos_zenith;
   Double_t        gandalf_dir_x;
   Double_t        gandalf_dir_y;
   Double_t        gandalf_dir_z;
   Double_t        gandalf_energy;
   Double_t        gandalf_id;
   Double_t        gandalf_is_good;
   Double_t        gandalf_length;
   Double_t        gandalf_likelihood;
   Double_t        gandalf_pos_x;
   Double_t        gandalf_pos_y;
   Double_t        gandalf_pos_z;
   Double_t        gandalf_rec_type;
   Double_t        gandalf_time;
   Double_t        shower_cos_zenith;
   Double_t        shower_dir_x;
   Double_t        shower_dir_y;
   Double_t        shower_dir_z;
   Double_t        shower_energy;
   Double_t        shower_id;
   Double_t        shower_is_good;
   Double_t        shower_length;
   Double_t        shower_likelihood;
   Double_t        shower_pos_x;
   Double_t        shower_pos_y;
   Double_t        shower_pos_z;
   Double_t        shower_rec_type;
   Double_t        shower_time;
   Double_t        dusj_DifferencesFirstAndSecondVertexFit_deltaTime;
   Double_t        dusj_DifferencesFirstAndSecondVertexFit_distance;
   Double_t        dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N;
   Double_t        dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference;
   Double_t        dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference;
   Double_t        dusj_Fork_muonSuppression_decision;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge2;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge_times_charge;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb2;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_charge;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_hit;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2__forXhitProb;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_hit2;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge;
   Double_t        dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_sum_XhitProb;
   Double_t        dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean;
   Double_t        dusj_L0AroundL1HitSelection__weight_withoutSelfSquared;
   Double_t        dusj_L0AroundL1HitSelection_weight;
   Double_t        dusj_L0AroundL1HitSelection_weight_withoutSelfSquared;
   Double_t        dusj_MuonSuppression_decision;
   Double_t        dusj_MuonSuppression_deltaTresQ20Q80;
   Double_t        dusj_MuonSuppression_enoughHits;
   Double_t        dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80;
   Double_t        dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75;
   Double_t        dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80;
   Double_t        dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Nom;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt_maxDeltaT10ns;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_azimuth;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorDown_bestLLH;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorUp_bestLLH;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy_bestLLH;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM__sum_forNoSignal;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM_sum;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum_forNoSignal;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llh_sum;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_multiplicity;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpFromPoissonOMhits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpOMhits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredOMhits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredPMThits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_zenith;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_azimuth;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_energy;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE0p8;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE1p2;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_overAllNorm;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum_forNoSignal;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_total;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_meanNomOverAllNorm;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_multiplicity;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_premiumEventFraction;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_pull;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_relativeWeightForOverAllNorm;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sigmaNomOverAllNorm;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpFromPoissonOMhits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpOMhits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredOMhits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredPMThits;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_zenith;
   Double_t        dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber;
   Double_t        dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber;
   Double_t        dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber;
   Double_t        dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber;
   Double_t        dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber;
   Double_t        dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N;
   Double_t        dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference;
   Double_t        dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N;
   Double_t        dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference;
   Double_t        dusj_dir_x;
   Double_t        dusj_dir_y;
   Double_t        dusj_dir_z;
   Double_t        dusj_energy;
   Double_t        dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;
   Double_t        dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;
   Double_t        dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;
   Double_t        dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;
   Double_t        dusj_id;
   Double_t        dusj_length;
   Double_t        dusj_likelihood;
   Double_t        dusj_pos_x;
   Double_t        dusj_pos_y;
   Double_t        dusj_pos_z;
   Double_t        dusj_rec_type;
   Double_t        dusj_time;
   Double_t        dusj_cos_zenith;
   Double_t        dusj_is_good;
   Double_t        group_id;
   Double_t        weight_w3;
   Double_t        weight_w2;
   Double_t        weight_w1;
   Double_t        run_id;
   Double_t        timestamp;
   Double_t        nanoseconds;
   Double_t        mc_time;
   Double_t        event_id;
   Double_t        mc_id;
   Double_t        is_neutrino;
   Double_t        dir_x;
   Double_t        dir_y;
   Double_t        dir_z;
   Double_t        pos_x;
   Double_t        pos_y;
   Double_t        pos_z;
   Double_t        energy;
   Double_t        bjorkeny;
   Double_t        interaction_channel;
   Double_t        is_cc;
   Double_t        type;
   Double_t        livetime_sec;
   Double_t        n_events_gen;
   Double_t        n_files_gen;
   Double_t        cos_zenith;
   Double_t        Erange_min;
   Double_t        Erange_max;
   Double_t        shower_upward_fraction;
   Double_t        shower_updown_deltaChi2;
   Double_t        gandalf_upward_fraction;
   Double_t        gandalf_updown_deltaChi2;
   Double_t        dusj_cosangle_dir_pos;
   Double_t        shower_cosangle_dir_pos;
   Double_t        gandalf_cosangle_dir_pos;
   Double_t        cosangle_dusj_shower;
   Double_t        cosangle_dusj_gandalf;
   Double_t        cosangle_shower_gandalf;
   Double_t        shower_shifted30m_pos_x;
   Double_t        shower_shifted30m_pos_y;
   Double_t        shower_shifted30m_pos_z;
   Double_t        gandalf_shifted30m_pos_x;
   Double_t        gandalf_shifted30m_pos_y;
   Double_t        gandalf_shifted30m_pos_z;
   Double_t        dusj_shifted30m_pos_x;
   Double_t        dusj_shifted30m_pos_y;
   Double_t        dusj_shifted30m_pos_z;
   Double_t        gandalf_lambda_antareslike;
   Double_t        weight_one_year;
   Double_t        _crossValidationTrainingMask_muon_score;
   Double_t        muon_score;
   Double_t        _crossValidationTrainingMask_track_score;
   Double_t        track_score;
   Double_t        _crossValidationTrainingMask_noise_score;
   Double_t        noise_score;
   Double_t        _crossValidationTrainingMask_muon_score_retrain;
   Double_t        muon_score_retrain;
   Double_t        dusj_energy_corrected;

   // List of branches
   TBranch        *b_gandalf_JENERGY_CHI2;   //!
   TBranch        *b_gandalf_JENERGY_ENERGY;   //!
   TBranch        *b_gandalf_JENERGY_MUON_RANGE_METRES;   //!
   TBranch        *b_gandalf_JENERGY_NDF;   //!
   TBranch        *b_gandalf_JENERGY_NOISE_LIKELIHOOD;   //!
   TBranch        *b_gandalf_JENERGY_NUMBER_OF_HITS;   //!
   TBranch        *b_gandalf_JGANDALF_BETA0_RAD;   //!
   TBranch        *b_gandalf_JGANDALF_BETA1_RAD;   //!
   TBranch        *b_gandalf_JGANDALF_CHI2;   //!
   TBranch        *b_gandalf_JGANDALF_LAMBDA;   //!
   TBranch        *b_gandalf_JGANDALF_NUMBER_OF_HITS;   //!
   TBranch        *b_gandalf_JGANDALF_NUMBER_OF_ITERATIONS;   //!
   TBranch        *b_gandalf_JSTART_LENGTH_METRES;   //!
   TBranch        *b_gandalf_JSTART_NPE_MIP;   //!
   TBranch        *b_gandalf_JSTART_NPE_MIP_TOTAL;   //!
   TBranch        *b_gandalf_JVETO_NPE;   //!
   TBranch        *b_gandalf_JVETO_NUMBER_OF_HITS;   //!
   TBranch        *b_gandalf_cos_zenith;   //!
   TBranch        *b_gandalf_dir_x;   //!
   TBranch        *b_gandalf_dir_y;   //!
   TBranch        *b_gandalf_dir_z;   //!
   TBranch        *b_gandalf_energy;   //!
   TBranch        *b_gandalf_id;   //!
   TBranch        *b_gandalf_is_good;   //!
   TBranch        *b_gandalf_length;   //!
   TBranch        *b_gandalf_likelihood;   //!
   TBranch        *b_gandalf_pos_x;   //!
   TBranch        *b_gandalf_pos_y;   //!
   TBranch        *b_gandalf_pos_z;   //!
   TBranch        *b_gandalf_rec_type;   //!
   TBranch        *b_gandalf_time;   //!
   TBranch        *b_shower_cos_zenith;   //!
   TBranch        *b_shower_dir_x;   //!
   TBranch        *b_shower_dir_y;   //!
   TBranch        *b_shower_dir_z;   //!
   TBranch        *b_shower_energy;   //!
   TBranch        *b_shower_id;   //!
   TBranch        *b_shower_is_good;   //!
   TBranch        *b_shower_length;   //!
   TBranch        *b_shower_likelihood;   //!
   TBranch        *b_shower_pos_x;   //!
   TBranch        *b_shower_pos_y;   //!
   TBranch        *b_shower_pos_z;   //!
   TBranch        *b_shower_rec_type;   //!
   TBranch        *b_shower_time;   //!
   TBranch        *b_dusj_DifferencesFirstAndSecondVertexFit_deltaTime;   //!
   TBranch        *b_dusj_DifferencesFirstAndSecondVertexFit_distance;   //!
   TBranch        *b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N;   //!
   TBranch        *b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference;   //!
   TBranch        *b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference;   //!
   TBranch        *b_dusj_Fork_muonSuppression_decision;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge2;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge_times_charge;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb2;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_charge;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_hit;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2__forXhitProb;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_hit2;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge;   //!
   TBranch        *b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_sum_XhitProb;   //!
   TBranch        *b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean;   //!
   TBranch        *b_dusj_L0AroundL1HitSelection__weight_withoutSelfSquared;   //!
   TBranch        *b_dusj_L0AroundL1HitSelection_weight;   //!
   TBranch        *b_dusj_L0AroundL1HitSelection_weight_withoutSelfSquared;   //!
   TBranch        *b_dusj_MuonSuppression_decision;   //!
   TBranch        *b_dusj_MuonSuppression_deltaTresQ20Q80;   //!
   TBranch        *b_dusj_MuonSuppression_enoughHits;   //!
   TBranch        *b_dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80;   //!
   TBranch        *b_dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75;   //!
   TBranch        *b_dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80;   //!
   TBranch        *b_dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Nom;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt_maxDeltaT10ns;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_azimuth;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorDown_bestLLH;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorUp_bestLLH;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy_bestLLH;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM__sum_forNoSignal;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM_sum;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum_forNoSignal;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llh_sum;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_multiplicity;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpFromPoissonOMhits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpOMhits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredOMhits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredPMThits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_zenith;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_azimuth;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_energy;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE0p8;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE1p2;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_overAllNorm;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum_forNoSignal;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_total;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_meanNomOverAllNorm;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_multiplicity;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_premiumEventFraction;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_pull;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_relativeWeightForOverAllNorm;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sigmaNomOverAllNorm;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpFromPoissonOMhits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpOMhits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredOMhits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredPMThits;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_zenith;   //!
   TBranch        *b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber;   //!
   TBranch        *b_dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber;   //!
   TBranch        *b_dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber;   //!
   TBranch        *b_dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber;   //!
   TBranch        *b_dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber;   //!
   TBranch        *b_dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N;   //!
   TBranch        *b_dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference;   //!
   TBranch        *b_dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N;   //!
   TBranch        *b_dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference;   //!
   TBranch        *b_dusj_dir_x;   //!
   TBranch        *b_dusj_dir_y;   //!
   TBranch        *b_dusj_dir_z;   //!
   TBranch        *b_dusj_energy;   //!
   TBranch        *b_dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;   //!
   TBranch        *b_dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;   //!
   TBranch        *b_dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;   //!
   TBranch        *b_dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult;   //!
   TBranch        *b_dusj_id;   //!
   TBranch        *b_dusj_length;   //!
   TBranch        *b_dusj_likelihood;   //!
   TBranch        *b_dusj_pos_x;   //!
   TBranch        *b_dusj_pos_y;   //!
   TBranch        *b_dusj_pos_z;   //!
   TBranch        *b_dusj_rec_type;   //!
   TBranch        *b_dusj_time;   //!
   TBranch        *b_dusj_cos_zenith;   //!
   TBranch        *b_dusj_is_good;   //!
   TBranch        *b_group_id;   //!
   TBranch        *b_weight_w3;   //!
   TBranch        *b_weight_w2;   //!
   TBranch        *b_weight_w1;   //!
   TBranch        *b_run_id;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_nanoseconds;   //!
   TBranch        *b_mc_time;   //!
   TBranch        *b_event_id;   //!
   TBranch        *b_mc_id;   //!
   TBranch        *b_is_neutrino;   //!
   TBranch        *b_dir_x;   //!
   TBranch        *b_dir_y;   //!
   TBranch        *b_dir_z;   //!
   TBranch        *b_pos_x;   //!
   TBranch        *b_pos_y;   //!
   TBranch        *b_pos_z;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_bjorkeny;   //!
   TBranch        *b_interaction_channel;   //!
   TBranch        *b_is_cc;   //!
   TBranch        *b_type;   //!
   TBranch        *b_livetime_sec;   //!
   TBranch        *b_n_events_gen;   //!
   TBranch        *b_n_files_gen;   //!
   TBranch        *b_cos_zenith;   //!
   TBranch        *b_Erange_min;   //!
   TBranch        *b_Erange_max;   //!
   TBranch        *b_shower_upward_fraction;   //!
   TBranch        *b_shower_updown_deltaChi2;   //!
   TBranch        *b_gandalf_upward_fraction;   //!
   TBranch        *b_gandalf_updown_deltaChi2;   //!
   TBranch        *b_dusj_cosangle_dir_pos;   //!
   TBranch        *b_shower_cosangle_dir_pos;   //!
   TBranch        *b_gandalf_cosangle_dir_pos;   //!
   TBranch        *b_cosangle_dusj_shower;   //!
   TBranch        *b_cosangle_dusj_gandalf;   //!
   TBranch        *b_cosangle_shower_gandalf;   //!
   TBranch        *b_shower_shifted30m_pos_x;   //!
   TBranch        *b_shower_shifted30m_pos_y;   //!
   TBranch        *b_shower_shifted30m_pos_z;   //!
   TBranch        *b_gandalf_shifted30m_pos_x;   //!
   TBranch        *b_gandalf_shifted30m_pos_y;   //!
   TBranch        *b_gandalf_shifted30m_pos_z;   //!
   TBranch        *b_dusj_shifted30m_pos_x;   //!
   TBranch        *b_dusj_shifted30m_pos_y;   //!
   TBranch        *b_dusj_shifted30m_pos_z;   //!
   TBranch        *b_gandalf_lambda_antareslike;   //!
   TBranch        *b_weight_one_year;   //!
   TBranch        *b__crossValidationTrainingMask_muon_score;   //!
   TBranch        *b_muon_score;   //!
   TBranch        *b__crossValidationTrainingMask_track_score;   //!
   TBranch        *b_track_score;   //!
   TBranch        *b__crossValidationTrainingMask_noise_score;   //!
   TBranch        *b_noise_score;   //!
   TBranch        *b__crossValidationTrainingMask_muon_score_retrain;   //!
   TBranch        *b_muon_score_retrain;   //!
   TBranch        *b_dusj_energy_corrected;   //!

   PIDBeta(TTree *tree=0);
   virtual ~PIDBeta();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PIDBeta_cxx
PIDBeta::PIDBeta(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../data/pid_result_13Oct2018_ORCA7.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../data/pid_result_13Oct2018_ORCA7.root");
      }
      f->GetObject("PID",tree);

   }
   Init(tree);
}

PIDBeta::~PIDBeta()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PIDBeta::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PIDBeta::LoadTree(Long64_t entry)
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

void PIDBeta::Init(TTree *tree)
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

   fChain->SetBranchAddress("gandalf_JENERGY_CHI2", &gandalf_JENERGY_CHI2, &b_gandalf_JENERGY_CHI2);
   fChain->SetBranchAddress("gandalf_JENERGY_ENERGY", &gandalf_JENERGY_ENERGY, &b_gandalf_JENERGY_ENERGY);
   fChain->SetBranchAddress("gandalf_JENERGY_MUON_RANGE_METRES", &gandalf_JENERGY_MUON_RANGE_METRES, &b_gandalf_JENERGY_MUON_RANGE_METRES);
   fChain->SetBranchAddress("gandalf_JENERGY_NDF", &gandalf_JENERGY_NDF, &b_gandalf_JENERGY_NDF);
   fChain->SetBranchAddress("gandalf_JENERGY_NOISE_LIKELIHOOD", &gandalf_JENERGY_NOISE_LIKELIHOOD, &b_gandalf_JENERGY_NOISE_LIKELIHOOD);
   fChain->SetBranchAddress("gandalf_JENERGY_NUMBER_OF_HITS", &gandalf_JENERGY_NUMBER_OF_HITS, &b_gandalf_JENERGY_NUMBER_OF_HITS);
   fChain->SetBranchAddress("gandalf_JGANDALF_BETA0_RAD", &gandalf_JGANDALF_BETA0_RAD, &b_gandalf_JGANDALF_BETA0_RAD);
   fChain->SetBranchAddress("gandalf_JGANDALF_BETA1_RAD", &gandalf_JGANDALF_BETA1_RAD, &b_gandalf_JGANDALF_BETA1_RAD);
   fChain->SetBranchAddress("gandalf_JGANDALF_CHI2", &gandalf_JGANDALF_CHI2, &b_gandalf_JGANDALF_CHI2);
   fChain->SetBranchAddress("gandalf_JGANDALF_LAMBDA", &gandalf_JGANDALF_LAMBDA, &b_gandalf_JGANDALF_LAMBDA);
   fChain->SetBranchAddress("gandalf_JGANDALF_NUMBER_OF_HITS", &gandalf_JGANDALF_NUMBER_OF_HITS, &b_gandalf_JGANDALF_NUMBER_OF_HITS);
   fChain->SetBranchAddress("gandalf_JGANDALF_NUMBER_OF_ITERATIONS", &gandalf_JGANDALF_NUMBER_OF_ITERATIONS, &b_gandalf_JGANDALF_NUMBER_OF_ITERATIONS);
   fChain->SetBranchAddress("gandalf_JSTART_LENGTH_METRES", &gandalf_JSTART_LENGTH_METRES, &b_gandalf_JSTART_LENGTH_METRES);
   fChain->SetBranchAddress("gandalf_JSTART_NPE_MIP", &gandalf_JSTART_NPE_MIP, &b_gandalf_JSTART_NPE_MIP);
   fChain->SetBranchAddress("gandalf_JSTART_NPE_MIP_TOTAL", &gandalf_JSTART_NPE_MIP_TOTAL, &b_gandalf_JSTART_NPE_MIP_TOTAL);
   fChain->SetBranchAddress("gandalf_JVETO_NPE", &gandalf_JVETO_NPE, &b_gandalf_JVETO_NPE);
   fChain->SetBranchAddress("gandalf_JVETO_NUMBER_OF_HITS", &gandalf_JVETO_NUMBER_OF_HITS, &b_gandalf_JVETO_NUMBER_OF_HITS);
   fChain->SetBranchAddress("gandalf_cos_zenith", &gandalf_cos_zenith, &b_gandalf_cos_zenith);
   fChain->SetBranchAddress("gandalf_dir_x", &gandalf_dir_x, &b_gandalf_dir_x);
   fChain->SetBranchAddress("gandalf_dir_y", &gandalf_dir_y, &b_gandalf_dir_y);
   fChain->SetBranchAddress("gandalf_dir_z", &gandalf_dir_z, &b_gandalf_dir_z);
   fChain->SetBranchAddress("gandalf_energy", &gandalf_energy, &b_gandalf_energy);
   fChain->SetBranchAddress("gandalf_id", &gandalf_id, &b_gandalf_id);
   fChain->SetBranchAddress("gandalf_is_good", &gandalf_is_good, &b_gandalf_is_good);
   fChain->SetBranchAddress("gandalf_length", &gandalf_length, &b_gandalf_length);
   fChain->SetBranchAddress("gandalf_likelihood", &gandalf_likelihood, &b_gandalf_likelihood);
   fChain->SetBranchAddress("gandalf_pos_x", &gandalf_pos_x, &b_gandalf_pos_x);
   fChain->SetBranchAddress("gandalf_pos_y", &gandalf_pos_y, &b_gandalf_pos_y);
   fChain->SetBranchAddress("gandalf_pos_z", &gandalf_pos_z, &b_gandalf_pos_z);
   fChain->SetBranchAddress("gandalf_rec_type", &gandalf_rec_type, &b_gandalf_rec_type);
   fChain->SetBranchAddress("gandalf_time", &gandalf_time, &b_gandalf_time);
   fChain->SetBranchAddress("shower_cos_zenith", &shower_cos_zenith, &b_shower_cos_zenith);
   fChain->SetBranchAddress("shower_dir_x", &shower_dir_x, &b_shower_dir_x);
   fChain->SetBranchAddress("shower_dir_y", &shower_dir_y, &b_shower_dir_y);
   fChain->SetBranchAddress("shower_dir_z", &shower_dir_z, &b_shower_dir_z);
   fChain->SetBranchAddress("shower_energy", &shower_energy, &b_shower_energy);
   fChain->SetBranchAddress("shower_id", &shower_id, &b_shower_id);
   fChain->SetBranchAddress("shower_is_good", &shower_is_good, &b_shower_is_good);
   fChain->SetBranchAddress("shower_length", &shower_length, &b_shower_length);
   fChain->SetBranchAddress("shower_likelihood", &shower_likelihood, &b_shower_likelihood);
   fChain->SetBranchAddress("shower_pos_x", &shower_pos_x, &b_shower_pos_x);
   fChain->SetBranchAddress("shower_pos_y", &shower_pos_y, &b_shower_pos_y);
   fChain->SetBranchAddress("shower_pos_z", &shower_pos_z, &b_shower_pos_z);
   fChain->SetBranchAddress("shower_rec_type", &shower_rec_type, &b_shower_rec_type);
   fChain->SetBranchAddress("shower_time", &shower_time, &b_shower_time);
   fChain->SetBranchAddress("dusj_DifferencesFirstAndSecondVertexFit_deltaTime", &dusj_DifferencesFirstAndSecondVertexFit_deltaTime, &b_dusj_DifferencesFirstAndSecondVertexFit_deltaTime);
   fChain->SetBranchAddress("dusj_DifferencesFirstAndSecondVertexFit_distance", &dusj_DifferencesFirstAndSecondVertexFit_distance, &b_dusj_DifferencesFirstAndSecondVertexFit_distance);
   fChain->SetBranchAddress("dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N", &dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N, &b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_N);
   fChain->SetBranchAddress("dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference", &dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference, &b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_difference);
   fChain->SetBranchAddress("dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference", &dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference, &b_dusj_FinalShowerHits_0dist60_L1cc_SingleHits_emisAng40_trackShowerTres_meanDifference);
   fChain->SetBranchAddress("dusj_Fork_muonSuppression_decision", &dusj_Fork_muonSuppression_decision, &b_dusj_Fork_muonSuppression_decision);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_correlationCoefficient);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_scalarProduct);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge2", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge2, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge2);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge_times_charge", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge_times_charge, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_Xcharge_times_charge);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb2", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb2, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb2);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_charge", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_charge, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_charge);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_hit", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_hit, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_XhitProb_times_hit);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2__forXhitProb", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2__forXhitProb, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_charge2__forXhitProb);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_hit2", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_hit2, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Charge_sum_hit2);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProbCharge_scalarProduct);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_HitProb_scalarProduct);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_NXcharge);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Ncharge);
   fChain->SetBranchAddress("dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits", &dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits, &b_dusj_HitPatternCharge_finalShowerHits_10degAroundCherAngle_Nhits);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_Mean);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_AbsTres_weightedMean);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Nhits);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_Mean);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_sum_XhitProb", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_sum_XhitProb, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_sum_XhitProb);
   fChain->SetBranchAddress("dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean", &dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean, &b_dusj_HitPatternTres_finalShowerHits_10degAroundCherAngle_Tres_weightedMean);
   fChain->SetBranchAddress("dusj_L0AroundL1HitSelection__weight_withoutSelfSquared", &dusj_L0AroundL1HitSelection__weight_withoutSelfSquared, &b_dusj_L0AroundL1HitSelection__weight_withoutSelfSquared);
   fChain->SetBranchAddress("dusj_L0AroundL1HitSelection_weight", &dusj_L0AroundL1HitSelection_weight, &b_dusj_L0AroundL1HitSelection_weight);
   fChain->SetBranchAddress("dusj_L0AroundL1HitSelection_weight_withoutSelfSquared", &dusj_L0AroundL1HitSelection_weight_withoutSelfSquared, &b_dusj_L0AroundL1HitSelection_weight_withoutSelfSquared);
   fChain->SetBranchAddress("dusj_MuonSuppression_decision", &dusj_MuonSuppression_decision, &b_dusj_MuonSuppression_decision);
   fChain->SetBranchAddress("dusj_MuonSuppression_deltaTresQ20Q80", &dusj_MuonSuppression_deltaTresQ20Q80, &b_dusj_MuonSuppression_deltaTresQ20Q80);
   fChain->SetBranchAddress("dusj_MuonSuppression_enoughHits", &dusj_MuonSuppression_enoughHits, &b_dusj_MuonSuppression_enoughHits);
   fChain->SetBranchAddress("dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80", &dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80, &b_dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80);
   fChain->SetBranchAddress("dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75", &dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75, &b_dusj_Trigger_3L1Dmax52_FinalShowerHits_0dist80m25tres75);
   fChain->SetBranchAddress("dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80", &dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80, &b_dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80);
   fChain->SetBranchAddress("dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75", &dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75, &b_dusj_Trigger_MX8hitsDmax46_FinalShowerHits_0dist80m25tres75);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Nom", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Nom, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Nom);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt_maxDeltaT10ns", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt_maxDeltaT10ns, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_Npmt_maxDeltaT10ns);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_azimuth", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_azimuth, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_azimuth);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorDown_bestLLH", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorDown_bestLLH, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorDown_bestLLH);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorUp_bestLLH", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorUp_bestLLH, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energyErrorUp_bestLLH);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy_bestLLH", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy_bestLLH, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_energy_bestLLH);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM__sum_forNoSignal", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM__sum_forNoSignal, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM__sum_forNoSignal);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM_sum", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM_sum, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhBestSinglePMTperDOM_sum);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum_forNoSignal", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum_forNoSignal, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llhSinglePMT_sum_forNoSignal);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llh_sum", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llh_sum, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_llh_sum);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_multiplicity", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_multiplicity, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_multiplicity);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpFromPoissonOMhits", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpFromPoissonOMhits, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpFromPoissonOMhits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpOMhits", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpOMhits, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumExpOMhits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredOMhits", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredOMhits, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredOMhits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredPMThits", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredPMThits, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_sumMeasuredPMThits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_zenith", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_zenith, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_FinalLLHValues_zenith);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_azimuth", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_azimuth, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_azimuth);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_energy", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_energy, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_energy);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE0p8", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE0p8, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE0p8);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE1p2", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE1p2, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_diff_forE1p2);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_overAllNorm", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_overAllNorm, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_overAllNorm);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum_forNoSignal", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum_forNoSignal, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_sum_forNoSignal);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_total", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_total, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_llh_total);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_meanNomOverAllNorm", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_meanNomOverAllNorm, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_meanNomOverAllNorm);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_multiplicity", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_multiplicity, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_multiplicity);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_premiumEventFraction", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_premiumEventFraction, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_premiumEventFraction);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_pull", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_pull, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_pull);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_relativeWeightForOverAllNorm", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_relativeWeightForOverAllNorm, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_relativeWeightForOverAllNorm);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sigmaNomOverAllNorm", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sigmaNomOverAllNorm, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sigmaNomOverAllNorm);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpFromPoissonOMhits", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpFromPoissonOMhits, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpFromPoissonOMhits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpOMhits", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpOMhits, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumExpOMhits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredOMhits", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredOMhits, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredOMhits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredPMThits", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredPMThits, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_sumMeasuredPMThits);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_zenith", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_zenith, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult_aroundCherAngle_FinalLLHValues_zenith);
   fChain->SetBranchAddress("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber", &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber, &b_dusj_best_DusjOrcaUsingProbabilitiesFinalFit_OUTVicinityNumber);
   fChain->SetBranchAddress("dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber", &dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber, &b_dusj_best_FirstDusjOrcaVertexFit_OUTVicinityNumber);
   fChain->SetBranchAddress("dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber", &dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber, &b_dusj_best_FirstDusjOrcaVertexFit_OUTVicinityWithTimeResidualToSeedNumber);
   fChain->SetBranchAddress("dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber", &dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber, &b_dusj_best_SecondDusjOrcaVertexFit_OUTFiducalNumber);
   fChain->SetBranchAddress("dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber", &dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber, &b_dusj_best_SecondDusjOrcaVertexFit_OUTVicinityNumber);
   fChain->SetBranchAddress("dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N", &dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N, &b_dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_N);
   fChain->SetBranchAddress("dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference", &dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference, &b_dusj_deltaTres_Q20_Q80_ClusteredL2ORV1L1HitSelection_SingleHits_difference);
   fChain->SetBranchAddress("dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N", &dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N, &b_dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_N);
   fChain->SetBranchAddress("dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference", &dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference, &b_dusj_deltaTres_Q20_Q80_FinalShowerHits_0dist60_L1cc_SingleHits_difference);
   fChain->SetBranchAddress("dusj_dir_x", &dusj_dir_x, &b_dusj_dir_x);
   fChain->SetBranchAddress("dusj_dir_y", &dusj_dir_y, &b_dusj_dir_y);
   fChain->SetBranchAddress("dusj_dir_z", &dusj_dir_z, &b_dusj_dir_z);
   fChain->SetBranchAddress("dusj_energy", &dusj_energy, &b_dusj_energy);
   fChain->SetBranchAddress("dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult", &dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult, &b_dusj_geoCoverage_R130h160_angle20_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult);
   fChain->SetBranchAddress("dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult", &dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult, &b_dusj_geoCoverage_R130h160_angle45_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult);
   fChain->SetBranchAddress("dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult", &dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult, &b_dusj_geoCoverage_R130h160_angle60_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult);
   fChain->SetBranchAddress("dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult", &dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult, &b_dusj_geoCoverage_R130h160_angle75_lmin30_best_DusjOrcaUsingProbabilitiesFinalFit_FitResult);
   fChain->SetBranchAddress("dusj_id", &dusj_id, &b_dusj_id);
   fChain->SetBranchAddress("dusj_length", &dusj_length, &b_dusj_length);
   fChain->SetBranchAddress("dusj_likelihood", &dusj_likelihood, &b_dusj_likelihood);
   fChain->SetBranchAddress("dusj_pos_x", &dusj_pos_x, &b_dusj_pos_x);
   fChain->SetBranchAddress("dusj_pos_y", &dusj_pos_y, &b_dusj_pos_y);
   fChain->SetBranchAddress("dusj_pos_z", &dusj_pos_z, &b_dusj_pos_z);
   fChain->SetBranchAddress("dusj_rec_type", &dusj_rec_type, &b_dusj_rec_type);
   fChain->SetBranchAddress("dusj_time", &dusj_time, &b_dusj_time);
   fChain->SetBranchAddress("dusj_cos_zenith", &dusj_cos_zenith, &b_dusj_cos_zenith);
   fChain->SetBranchAddress("dusj_is_good", &dusj_is_good, &b_dusj_is_good);
   fChain->SetBranchAddress("group_id", &group_id, &b_group_id);
   fChain->SetBranchAddress("weight_w3", &weight_w3, &b_weight_w3);
   fChain->SetBranchAddress("weight_w2", &weight_w2, &b_weight_w2);
   fChain->SetBranchAddress("weight_w1", &weight_w1, &b_weight_w1);
   fChain->SetBranchAddress("run_id", &run_id, &b_run_id);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("nanoseconds", &nanoseconds, &b_nanoseconds);
   fChain->SetBranchAddress("mc_time", &mc_time, &b_mc_time);
   fChain->SetBranchAddress("event_id", &event_id, &b_event_id);
   fChain->SetBranchAddress("mc_id", &mc_id, &b_mc_id);
   fChain->SetBranchAddress("is_neutrino", &is_neutrino, &b_is_neutrino);
   fChain->SetBranchAddress("dir_x", &dir_x, &b_dir_x);
   fChain->SetBranchAddress("dir_y", &dir_y, &b_dir_y);
   fChain->SetBranchAddress("dir_z", &dir_z, &b_dir_z);
   fChain->SetBranchAddress("pos_x", &pos_x, &b_pos_x);
   fChain->SetBranchAddress("pos_y", &pos_y, &b_pos_y);
   fChain->SetBranchAddress("pos_z", &pos_z, &b_pos_z);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("bjorkeny", &bjorkeny, &b_bjorkeny);
   fChain->SetBranchAddress("interaction_channel", &interaction_channel, &b_interaction_channel);
   fChain->SetBranchAddress("is_cc", &is_cc, &b_is_cc);
   fChain->SetBranchAddress("type", &type, &b_type);
   fChain->SetBranchAddress("livetime_sec", &livetime_sec, &b_livetime_sec);
   fChain->SetBranchAddress("n_events_gen", &n_events_gen, &b_n_events_gen);
   fChain->SetBranchAddress("n_files_gen", &n_files_gen, &b_n_files_gen);
   fChain->SetBranchAddress("cos_zenith", &cos_zenith, &b_cos_zenith);
   fChain->SetBranchAddress("Erange_min", &Erange_min, &b_Erange_min);
   fChain->SetBranchAddress("Erange_max", &Erange_max, &b_Erange_max);
   fChain->SetBranchAddress("shower_upward_fraction", &shower_upward_fraction, &b_shower_upward_fraction);
   fChain->SetBranchAddress("shower_updown_deltaChi2", &shower_updown_deltaChi2, &b_shower_updown_deltaChi2);
   fChain->SetBranchAddress("gandalf_upward_fraction", &gandalf_upward_fraction, &b_gandalf_upward_fraction);
   fChain->SetBranchAddress("gandalf_updown_deltaChi2", &gandalf_updown_deltaChi2, &b_gandalf_updown_deltaChi2);
   fChain->SetBranchAddress("dusj_cosangle_dir_pos", &dusj_cosangle_dir_pos, &b_dusj_cosangle_dir_pos);
   fChain->SetBranchAddress("shower_cosangle_dir_pos", &shower_cosangle_dir_pos, &b_shower_cosangle_dir_pos);
   fChain->SetBranchAddress("gandalf_cosangle_dir_pos", &gandalf_cosangle_dir_pos, &b_gandalf_cosangle_dir_pos);
   fChain->SetBranchAddress("cosangle_dusj_shower", &cosangle_dusj_shower, &b_cosangle_dusj_shower);
   fChain->SetBranchAddress("cosangle_dusj_gandalf", &cosangle_dusj_gandalf, &b_cosangle_dusj_gandalf);
   fChain->SetBranchAddress("cosangle_shower_gandalf", &cosangle_shower_gandalf, &b_cosangle_shower_gandalf);
   fChain->SetBranchAddress("shower_shifted30m_pos_x", &shower_shifted30m_pos_x, &b_shower_shifted30m_pos_x);
   fChain->SetBranchAddress("shower_shifted30m_pos_y", &shower_shifted30m_pos_y, &b_shower_shifted30m_pos_y);
   fChain->SetBranchAddress("shower_shifted30m_pos_z", &shower_shifted30m_pos_z, &b_shower_shifted30m_pos_z);
   fChain->SetBranchAddress("gandalf_shifted30m_pos_x", &gandalf_shifted30m_pos_x, &b_gandalf_shifted30m_pos_x);
   fChain->SetBranchAddress("gandalf_shifted30m_pos_y", &gandalf_shifted30m_pos_y, &b_gandalf_shifted30m_pos_y);
   fChain->SetBranchAddress("gandalf_shifted30m_pos_z", &gandalf_shifted30m_pos_z, &b_gandalf_shifted30m_pos_z);
   fChain->SetBranchAddress("dusj_shifted30m_pos_x", &dusj_shifted30m_pos_x, &b_dusj_shifted30m_pos_x);
   fChain->SetBranchAddress("dusj_shifted30m_pos_y", &dusj_shifted30m_pos_y, &b_dusj_shifted30m_pos_y);
   fChain->SetBranchAddress("dusj_shifted30m_pos_z", &dusj_shifted30m_pos_z, &b_dusj_shifted30m_pos_z);
   fChain->SetBranchAddress("gandalf_lambda_antareslike", &gandalf_lambda_antareslike, &b_gandalf_lambda_antareslike);
   fChain->SetBranchAddress("weight_one_year", &weight_one_year, &b_weight_one_year);
   fChain->SetBranchAddress("_crossValidationTrainingMask_muon_score", &_crossValidationTrainingMask_muon_score, &b__crossValidationTrainingMask_muon_score);
   fChain->SetBranchAddress("muon_score", &muon_score, &b_muon_score);
   fChain->SetBranchAddress("_crossValidationTrainingMask_track_score", &_crossValidationTrainingMask_track_score, &b__crossValidationTrainingMask_track_score);
   fChain->SetBranchAddress("track_score", &track_score, &b_track_score);
   fChain->SetBranchAddress("_crossValidationTrainingMask_noise_score", &_crossValidationTrainingMask_noise_score, &b__crossValidationTrainingMask_noise_score);
   fChain->SetBranchAddress("noise_score", &noise_score, &b_noise_score);
   fChain->SetBranchAddress("_crossValidationTrainingMask_muon_score_retrain", &_crossValidationTrainingMask_muon_score_retrain, &b__crossValidationTrainingMask_muon_score_retrain);
   fChain->SetBranchAddress("muon_score_retrain", &muon_score_retrain, &b_muon_score_retrain);
   fChain->SetBranchAddress("dusj_energy_corrected", &dusj_energy_corrected, &b_dusj_energy_corrected);
   Notify();
}

Bool_t PIDBeta::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PIDBeta::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PIDBeta::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PIDBeta_cxx
