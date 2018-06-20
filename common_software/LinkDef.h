#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//when you do root, .L libnmhsoft.so, you get access to the classes listed below.
//If you run a compiled macro and include before it does not matter
#pragma link C++ enum AtmFluxOpt;
#pragma link C++ class GSGParser+;
#pragma link C++ class SummaryParser+;
#pragma link C++ class AtmFlux+;
#pragma link C++ class NuXsec+;
#pragma link C++ class FileHeader+;
#pragma link C++ class SummaryEvent+;
#pragma link C++ class EventFilter+;
#pragma link C++ class EventSelection+;


#endif
