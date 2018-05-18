#include "GSGSampler.C"
using namespace std;

void evtid_functionality() {

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");

  //------------------------------------------------------
  // generate some pseudo-events
  //------------------------------------------------------

  TRandom3 rand(10);
  vector<evtid> evts;

  for (Int_t i = 0; i < 100; i++) {
    Int_t run_nr = rand.Uniform(0,10);
    Int_t evt_nr = rand.Uniform(0,100);
    Int_t e_min = 1;
    if ( rand.Uniform(0,1) > 0.5 ) e_min = 3;
    evts.push_back( evtid(run_nr, evt_nr, e_min) );
  }

  cout << "******************************************************" << endl;
  cout << "Before sort" << endl;
  cout << "******************************************************" << endl;

  Int_t idx = 0;
  for (auto evt: evts) {
    cout << idx << "\t" << evt.run_nr << "\t" << evt.evt_nr << "\t" << evt.e_min << endl;
    idx++;
  }

  cout << "******************************************************" << endl;
  cout << "After sort" << endl;
  cout << "******************************************************" << endl;

  std::sort( evts.begin(), evts.end() );

  idx = 0;
  for (auto evt: evts) {
    cout << idx << "\t" << evt.run_nr << "\t" << evt.evt_nr << "\t" << evt.e_min << endl;
  }

  cout << "******************************************************" << endl;
  cout << "Find" << endl;
  cout << "******************************************************" << endl;

  std::pair<Int_t, Int_t> range = get_start_stop( evts, 5, 3 );

  for (auto it = evts.begin() + range.first; it < evts.begin() + range.second; it++) {
    cout << it->run_nr << "\t" << it->evt_nr << "\t" << it->e_min << endl;
  }

  for (Int_t i = range.first; i < range.second; i++) {
    cout << i << "\t" << evts[i].run_nr << "\t" << evts[i].evt_nr << "\t" << evts[i].e_min << endl;
  }

  evtid test(6, 2, 3);
  if ( std::find(evts.begin() + range.first, evts.begin() + range.second, test) != ( evts.begin() + range.second) ) {
    cout << "Found" << endl;
  }
  else {
    cout << "Not found" << endl; 
  }

}
