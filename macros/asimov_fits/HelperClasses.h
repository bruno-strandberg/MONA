#ifndef HelperClasses_h
#define HelperClasses_h

#include "HelperFunctions.h"

#include "DetResponse.h"
#include "EventSelection.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <map>
#include <vector>
#include <tuple>

class DetResponseCreator {

 public:
  DetResponseCreator(Int_t n_pid_cats, TString path_to_file) {
      
    // Set default values for PID boundaries at given category.
    fNPidCategories = n_pid_cats;
    fPidMap = SetPIDCase(fNPidCategories);
    fDetResFilePath = SetFilePath(path_to_file);

    // Add the DRs to the vector for the PID category
    std::tuple< std::vector<DetResponse*>, std::vector<DetResponse*> > drWCuts = ResponseCreator();
    fTrkRespVector = std::get<1>(drWCuts);
    fShwRespVector = std::get<0>(drWCuts);

    fDetRespFilesExist = ResponseFilesExist();
  };

  DetResponseCreator(Int_t n_pid_cats, TString path_to_file, 
          Int_t ebins,  Double_t emin,  Double_t emax, 
          Int_t ctbins, Double_t ctmin, Double_t ctmax, 
          Int_t bybins, Double_t bymin, Double_t bymax) {

    // Set default values for PID boundaries at given category.
    fNPidCategories = n_pid_cats;
    fPidMap = SetPIDCase(fNPidCategories);
    fDetResFilePath = SetFilePath(path_to_file);

    // Overwrite the default binnings of the ResponseCreator
    std::tuple< std::vector<DetResponse*>, std::vector<DetResponse*> > drWCuts = ResponseCreator(ebins, emin, emax,
            ctbins, ctmin, ctmax, bybins, bymin, bymax);
    fTrkRespVector = std::get<1>(drWCuts);
    fShwRespVector = std::get<0>(drWCuts);

    fDetRespFilesExist = ResponseFilesExist();
  };

  ~DetResponseCreator() {
    for (Int_t i = 0; i < fNPidCategories; i++) {
      delete fTrkRespVector[i];
      delete fShwRespVector[i];
    }
  };


  TString GetFilePath() { return fDetResFilePath; }

  Bool_t GetFilesExist() { 
    if (fDetRespFilesExist) {
      cout << "NOTICE: The DR output files are available at " << fDetResFilePath << endl; 
    }
    else {
      cout << "NOTICE: The DR output files are not available" << endl; 
    }
    return fDetRespFilesExist;
  }

  void FillResponseVectors(SummaryParser* sumparser) {
    Bool_t files_exist = GetFilesExist();

    if (not files_exist) {
      SumParserFillResponseVectors(sumparser);
      WriteDetResFiles();
    }
    else {
      ReadDetResFiles();
    }
  }

  std::vector<DetResponse*> GetTrkResponseVector() { return fTrkRespVector; }
  std::vector<DetResponse*> GetShwResponseVector() { return fShwRespVector; }


 private:
  std::tuple< std::vector<DetResponse*>, std::vector<DetResponse*> > ResponseCreator( 
          Int_t ebins  = 24, Double_t emin  =  2., Double_t emax  = 80., 
          Int_t ctbins = 40, Double_t ctmin = -1., Double_t ctmax = 0., 
          Int_t bybins = 1 , Double_t bymin =  0., Double_t bymax = 1.) {
      
    for (Int_t i = 0; i < fNPidCategories; i++) {
      std::function<bool(double, double)> comparison_operator;
      if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first category needs to include the lower limit.
      } else { comparison_operator = std::greater<double>(); }

      DetResponse* fTrkResp = new DetResponse(DetResponse::track, Form("track_response_%.2f", fPidMap[i]),
                                 ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax);
      fTrkResp->AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5         , true );
      fTrkResp->AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5         , true );
      fTrkResp->AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , fPidMap[i]  , true );
      fTrkResp->AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), fPidMap[i+1], true );
      fTrkResp->AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05        , true );
      fTrkResp->AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18        , true );
      fTrkRespVector.push_back(fTrkResp);

      DetResponse* fShwResp = new DetResponse(DetResponse::shower, Form("shower_response_%.2f", fPidMap[i]),
                                  ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax);
      fShwResp->AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5         , true );
      fShwResp->AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5         , true );
      fShwResp->AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , fPidMap[i]  , true );
      fShwResp->AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), fPidMap[i+1], true );
      fShwResp->AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05        , true );
      fShwResp->AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5         , true );
      fShwRespVector.push_back(fShwResp);
    }

    return std::make_tuple(fShwRespVector, fTrkRespVector);
  };

  Bool_t ResponseFilesExist() {
    for (Int_t i = 0; i < fNPidCategories; i++) {
      TString track_file  = Form("track_response_%.2f.root" , fPidMap[i]);
      TString shower_file = Form("shower_response_%.2f.root", fPidMap[i]);
      Bool_t track_exists  = NMHUtils::FileExists(fDetResFilePath + track_file);
      Bool_t shower_exists = NMHUtils::FileExists(fDetResFilePath + shower_file);
      fDetRespFilesExist = ((kTRUE and track_exists) and shower_exists);
    }
    return fDetRespFilesExist;
  };

  void SumParserFillResponseVectors(SummaryParser* sumparser) {

    for (Int_t i = 0; i < sumparser->GetTree()->GetEntries(); i++) {
      if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
      SummaryEvent* evt = sumparser->GetEvt(i);
      for (Int_t j = 0; j < fNPidCategories; j++) {
        fTrkRespVector[j]->Fill(evt);
        fShwRespVector[j]->Fill(evt);
      }
    }
    cout << "NOTICE: Finished filling response through SummaryParser" << endl;
  }

  void WriteDetResFiles(TString path_to_file = "") {
    cout << "NOTICE: Writing responses to disk" << endl;

    if (path_to_file == "") { path_to_file = fDetResFilePath; }

    for (Int_t i = 0; i < fNPidCategories; i++) {
      fTrkRespVector[i]->WriteToFile( path_to_file + Form("track_response_%.2f.root",  fPidMap[i]) );
      fShwRespVector[i]->WriteToFile( path_to_file + Form("shower_response_%.2f.root", fPidMap[i]) );
    }
  }

  void ReadDetResFiles(TString path_to_file = "") {
    cout << "NOTICE: Reading responses from disk" << endl;
    
    if (path_to_file == "") { path_to_file = fDetResFilePath; }

    for (Int_t i = 0; i < fNPidCategories; i++) {
      fTrkRespVector[i]->ReadFromFile( path_to_file + Form("track_response_%.2f.root" , fPidMap[i]) );
      fShwRespVector[i]->ReadFromFile( path_to_file + Form("shower_response_%.2f.root", fPidMap[i]) );
    }
    cout << "NOTICE: Finished filling response through reading from disk" << endl;
  }

  TString SetFilePath(TString path_to_file) {
    TString filepath = path_to_file;
    Int_t stringLength = filepath.Length();
    if ( filepath[stringLength-1] != (TString)"/" ) filepath.Append("/");
    return filepath;
  }  

  TString fDetResFilePath;                  // File path where the DR root files are stored.
  Int_t fNPidCategories;                    // Number of PID categories.
  std::map<Int_t, Double_t> fPidMap;        // Map of the PID boundaries in Q-space, track-like quality space.
  std::vector<DetResponse*> fTrkRespVector; // Vector containing pointers to the the track responses.
  std::vector<DetResponse*> fShwRespVector; // Vector containing pointers to the the shower responses.
  Bool_t fDetRespFilesExist = kFALSE;       // Bool for the existence of the input/output DR root files.

};

#endif
