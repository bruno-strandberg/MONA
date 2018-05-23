#include "../evt_sampler/GSGSampler.C"
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
    cout << "Index: " << idx << " Event: " << evt << endl;
    idx++;
  }

  cout << "******************************************************" << endl;
  cout << "After sort" << endl;
  cout << "******************************************************" << endl;

  std::sort( evts.begin(), evts.end() );

  idx = 0;
  for (auto evt: evts) {
    cout << "Index: " << idx << " Event: " << evt << endl;
    idx++;
  }

  cout << "******************************************************" << endl;
  cout << "Search range" << endl;
  cout << "******************************************************" << endl;

  evtid test(5, 2, 3);
  std::pair<Int_t, Int_t> range = get_start_stop( evts, test.run_nr, test.e_min );

  for (auto it = evts.begin() + range.first; it < evts.begin() + range.second; it++) {
    cout << it->run_nr << "\t" << it->evt_nr << "\t" << it->e_min << endl;
  }

  for (Int_t i = range.first; i < range.second; i++) {
    cout << i << "\t" << evts[i].run_nr << "\t" << evts[i].evt_nr << "\t" << evts[i].e_min << endl;
  }

  cout << "******************************************************" << endl;
  cout << "Find" << endl;
  cout << "******************************************************" << endl;

  TStopwatch find_time;
  if ( std::find(evts.begin() + range.first, evts.begin() + range.second, test) != ( evts.begin() + range.second) ) {
    cout << "std:find: found" << endl;
  }
  else {
    cout << "std:find: not found" << endl;
  }
  cout << "std::find() time" << find_time.RealTime() << endl;
  

  TStopwatch bs_time;
  if ( std::binary_search(evts.begin() + range.first, evts.begin() + range.second, test) ) {
    cout << "Binary search: found" << endl;
  }
  else {
    cout << "Binary search: not found" << endl; 
  }
  cout << "binary search time" << bs_time.RealTime() << endl;

}
