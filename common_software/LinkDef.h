#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//when you do root, .L libcommonsoft.so, you get access to GSGParser and SummaryParser
//if you run a compiled macro and include before it does not matter
#pragma link C++ class GSGParser+;
#pragma link C++ class SummaryParser+;

#endif
