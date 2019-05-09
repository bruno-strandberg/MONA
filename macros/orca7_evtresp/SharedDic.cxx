// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME SharedDic

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "ORCA7.h"

// Header files passed via #pragma extra_include

namespace O7 {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *O7_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("O7", 0 /*version*/, "ORCA7.h", 20,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &O7_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *O7_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static TClass *ORCA7_Dictionary();
   static void ORCA7_TClassManip(TClass*);
   static void delete_ORCA7(void *p);
   static void deleteArray_ORCA7(void *p);
   static void destruct_ORCA7(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ORCA7*)
   {
      ::ORCA7 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ORCA7));
      static ::ROOT::TGenericClassInfo 
         instance("ORCA7", "ORCA7.h", 53,
                  typeid(::ORCA7), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ORCA7_Dictionary, isa_proxy, 4,
                  sizeof(::ORCA7) );
      instance.SetDelete(&delete_ORCA7);
      instance.SetDeleteArray(&deleteArray_ORCA7);
      instance.SetDestructor(&destruct_ORCA7);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ORCA7*)
   {
      return GenerateInitInstanceLocal((::ORCA7*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ORCA7*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ORCA7_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ORCA7*)0x0)->GetClass();
      ORCA7_TClassManip(theClass);
   return theClass;
   }

   static void ORCA7_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ORCA7(void *p) {
      delete ((::ORCA7*)p);
   }
   static void deleteArray_ORCA7(void *p) {
      delete [] ((::ORCA7*)p);
   }
   static void destruct_ORCA7(void *p) {
      typedef ::ORCA7 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ORCA7

namespace ROOT {
   static TClass *vectorlEO7cLcLPidBinConfgR_Dictionary();
   static void vectorlEO7cLcLPidBinConfgR_TClassManip(TClass*);
   static void *new_vectorlEO7cLcLPidBinConfgR(void *p = 0);
   static void *newArray_vectorlEO7cLcLPidBinConfgR(Long_t size, void *p);
   static void delete_vectorlEO7cLcLPidBinConfgR(void *p);
   static void deleteArray_vectorlEO7cLcLPidBinConfgR(void *p);
   static void destruct_vectorlEO7cLcLPidBinConfgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<O7::PidBinConf>*)
   {
      vector<O7::PidBinConf> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<O7::PidBinConf>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<O7::PidBinConf>", -2, "vector", 210,
                  typeid(vector<O7::PidBinConf>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEO7cLcLPidBinConfgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<O7::PidBinConf>) );
      instance.SetNew(&new_vectorlEO7cLcLPidBinConfgR);
      instance.SetNewArray(&newArray_vectorlEO7cLcLPidBinConfgR);
      instance.SetDelete(&delete_vectorlEO7cLcLPidBinConfgR);
      instance.SetDeleteArray(&deleteArray_vectorlEO7cLcLPidBinConfgR);
      instance.SetDestructor(&destruct_vectorlEO7cLcLPidBinConfgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<O7::PidBinConf> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<O7::PidBinConf>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEO7cLcLPidBinConfgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<O7::PidBinConf>*)0x0)->GetClass();
      vectorlEO7cLcLPidBinConfgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEO7cLcLPidBinConfgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEO7cLcLPidBinConfgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<O7::PidBinConf> : new vector<O7::PidBinConf>;
   }
   static void *newArray_vectorlEO7cLcLPidBinConfgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<O7::PidBinConf>[nElements] : new vector<O7::PidBinConf>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEO7cLcLPidBinConfgR(void *p) {
      delete ((vector<O7::PidBinConf>*)p);
   }
   static void deleteArray_vectorlEO7cLcLPidBinConfgR(void *p) {
      delete [] ((vector<O7::PidBinConf>*)p);
   }
   static void destruct_vectorlEO7cLcLPidBinConfgR(void *p) {
      typedef vector<O7::PidBinConf> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<O7::PidBinConf>

namespace ROOT {
   static TClass *maplETStringcOFitPDFmUgR_Dictionary();
   static void maplETStringcOFitPDFmUgR_TClassManip(TClass*);
   static void *new_maplETStringcOFitPDFmUgR(void *p = 0);
   static void *newArray_maplETStringcOFitPDFmUgR(Long_t size, void *p);
   static void delete_maplETStringcOFitPDFmUgR(void *p);
   static void deleteArray_maplETStringcOFitPDFmUgR(void *p);
   static void destruct_maplETStringcOFitPDFmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<TString,FitPDF*>*)
   {
      map<TString,FitPDF*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<TString,FitPDF*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<TString,FitPDF*>", -2, "map", 96,
                  typeid(map<TString,FitPDF*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplETStringcOFitPDFmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<TString,FitPDF*>) );
      instance.SetNew(&new_maplETStringcOFitPDFmUgR);
      instance.SetNewArray(&newArray_maplETStringcOFitPDFmUgR);
      instance.SetDelete(&delete_maplETStringcOFitPDFmUgR);
      instance.SetDeleteArray(&deleteArray_maplETStringcOFitPDFmUgR);
      instance.SetDestructor(&destruct_maplETStringcOFitPDFmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<TString,FitPDF*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<TString,FitPDF*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplETStringcOFitPDFmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<TString,FitPDF*>*)0x0)->GetClass();
      maplETStringcOFitPDFmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplETStringcOFitPDFmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplETStringcOFitPDFmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<TString,FitPDF*> : new map<TString,FitPDF*>;
   }
   static void *newArray_maplETStringcOFitPDFmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<TString,FitPDF*>[nElements] : new map<TString,FitPDF*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplETStringcOFitPDFmUgR(void *p) {
      delete ((map<TString,FitPDF*>*)p);
   }
   static void deleteArray_maplETStringcOFitPDFmUgR(void *p) {
      delete [] ((map<TString,FitPDF*>*)p);
   }
   static void destruct_maplETStringcOFitPDFmUgR(void *p) {
      typedef map<TString,FitPDF*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<TString,FitPDF*>

namespace ROOT {
   static TClass *maplETStringcOEvtResponsemUgR_Dictionary();
   static void maplETStringcOEvtResponsemUgR_TClassManip(TClass*);
   static void *new_maplETStringcOEvtResponsemUgR(void *p = 0);
   static void *newArray_maplETStringcOEvtResponsemUgR(Long_t size, void *p);
   static void delete_maplETStringcOEvtResponsemUgR(void *p);
   static void deleteArray_maplETStringcOEvtResponsemUgR(void *p);
   static void destruct_maplETStringcOEvtResponsemUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<TString,EvtResponse*>*)
   {
      map<TString,EvtResponse*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<TString,EvtResponse*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<TString,EvtResponse*>", -2, "map", 96,
                  typeid(map<TString,EvtResponse*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplETStringcOEvtResponsemUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<TString,EvtResponse*>) );
      instance.SetNew(&new_maplETStringcOEvtResponsemUgR);
      instance.SetNewArray(&newArray_maplETStringcOEvtResponsemUgR);
      instance.SetDelete(&delete_maplETStringcOEvtResponsemUgR);
      instance.SetDeleteArray(&deleteArray_maplETStringcOEvtResponsemUgR);
      instance.SetDestructor(&destruct_maplETStringcOEvtResponsemUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<TString,EvtResponse*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<TString,EvtResponse*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplETStringcOEvtResponsemUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<TString,EvtResponse*>*)0x0)->GetClass();
      maplETStringcOEvtResponsemUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplETStringcOEvtResponsemUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplETStringcOEvtResponsemUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<TString,EvtResponse*> : new map<TString,EvtResponse*>;
   }
   static void *newArray_maplETStringcOEvtResponsemUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<TString,EvtResponse*>[nElements] : new map<TString,EvtResponse*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplETStringcOEvtResponsemUgR(void *p) {
      delete ((map<TString,EvtResponse*>*)p);
   }
   static void deleteArray_maplETStringcOEvtResponsemUgR(void *p) {
      delete [] ((map<TString,EvtResponse*>*)p);
   }
   static void destruct_maplETStringcOEvtResponsemUgR(void *p) {
      typedef map<TString,EvtResponse*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<TString,EvtResponse*>

namespace {
  void TriggerDictionaryInitialization_SharedDic_Impl() {
    static const char* headers[] = {
"ORCA7.h",
0
    };
    static const char* includePaths[] = {
"/pbs/software/centos-7-x86_64/root/6.10.02/include/root",
"/sps/km3net/users/jmanczak/MONA/resp_ebe_new/macros/orca7/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "SharedDic dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate("$clingAutoload$ORCA7.h")))  ORCA7;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "SharedDic dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "ORCA7.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"ORCA7", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("SharedDic",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_SharedDic_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_SharedDic_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_SharedDic() {
  TriggerDictionaryInitialization_SharedDic_Impl();
}
