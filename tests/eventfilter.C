void eventfilter() {

  //-----------------------------------------------------------
  //create an event and add some random values for filtering
  //-----------------------------------------------------------

  SummaryEvent *evt = new SummaryEvent();
  evt->Set_track_energy(10);
  evt->Set_track_dir(0  , 0, -0.7);
  evt->Set_track_pos(-10, 0, -90);
  evt->Set_track_ql0(1);
  evt->Set_track_ql1(1);
  evt->Set_track_ql2(0);
  evt->Set_RDF_track_score(0.8);
  evt->Set_RDF_muon_score(0.1);
  evt->Set_RDF_noise_score(0.02);

  //-----------------------------------------------------------
  //create first filter, event expects > 5 GeV energy and up-going direction
  //-----------------------------------------------------------
  EventFilter f1;
  f1.AddCut(&SummaryEvent::Get_track_energy, 5, std::greater<double>(), true );
  f1.AddCut(&SummaryEvent::Get_track_dir_z, 0, std::less<double>(), true );

  cout << "Filter f1 (expect pass) " << f1.PassesCuts(evt) << endl;

  //-----------------------------------------------------------
  //create second filter, additionally as x or y direction > 0
  //-----------------------------------------------------------
  EventFilter f2;
  f2.AddCut(&SummaryEvent::Get_track_energy, 5, std::greater<double>(), true );
  f2.AddCut(&SummaryEvent::Get_track_dir_z, 0, std::less<double>(), true );

  f2.AddCut(&SummaryEvent::Get_track_dir_x, 0, std::greater<double>(), false );
  f2.AddCut(&SummaryEvent::Get_track_dir_y, 0, std::greater<double>(), false );
  
  cout << "Filter f2 (expect fail) " << f2.PassesCuts(evt) << endl;

}
