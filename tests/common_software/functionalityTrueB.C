#include "DetResponse.h"

#include <iostream>
using namespace std;

/** This application tests the functionality of `TrueB` structure, defined in `DetResponse.h`
    \return 0 if test worked, 1 otherwise.
*/

int main(const int argc, const char** argv) {

  // default constructor; check that everything is 0
  //-------------------------------------------------
  TrueB tb;

  if (tb.fFlav != 0 || tb.fIsCC != 0 || tb.fIsNB != 0 || tb.fE_true_bin != 0 || tb.fCt_true_bin != 0 || tb.fBy_true_bin != 0 || tb.fW != 0 || tb.fWE != 0) {
    cout << "NOTICE functionalityTrueB default constructor fails" << endl;
    return 1;
  }

  // try to construct with wrong flavor
  //-------------------------------------------------
  bool caught_err = false;

  try {
    TrueB(3,1,2,0,0,0);
  }
  catch (const std::invalid_argument& ia) {
    caught_err = true;
  }
  
  if ( !caught_err ) return 1;

  // test copy constructor
  //-------------------------------------------------
  TrueB O(2,1,0, 32, 12, 1);
  TrueB C( O );

  if (O.fFlav != C.fFlav || O.fIsCC != C.fIsCC || O.fIsNB != C.fIsNB || O.fE_true_bin != C.fE_true_bin ||
      O.fCt_true_bin != C.fCt_true_bin || O.fBy_true_bin != C.fBy_true_bin || O.fW != C.fW || O.fWE != C.fWE) {
    cout << "NOTICE functionalityTrueB copy constructor fails" << endl;
    return 1;
  }

  // test comparison operator
  //-------------------------------------------------
  if ( O != C ) {
    cout << "NOTICE functionalityTrueB comparison operator fails" << endl;
    return 1;
  }

  TrueB CC( C );
  CC.fFlav = 1;

  if ( C == CC ) {
    cout << "NOTICE functionalityTrueB comparison operator fails" << endl;
    return 1;
  }
  
  return 0;
  
}
