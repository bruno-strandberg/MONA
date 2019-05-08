#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//when you do root, .L libnmhsoft.so, you get access to the classes listed below.
//If you run a compiled macro and include before it does not matter
#pragma link C++ namespace O7+;
#pragma link C++ class fitpacket+;
#pragma link C++ class ORCA7+;

#endif
