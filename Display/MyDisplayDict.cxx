// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MyDisplayDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "MyDisplay.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_MyDisplay(void *p = nullptr);
   static void *newArray_MyDisplay(Long_t size, void *p);
   static void delete_MyDisplay(void *p);
   static void deleteArray_MyDisplay(void *p);
   static void destruct_MyDisplay(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MyDisplay*)
   {
      ::MyDisplay *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MyDisplay >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("MyDisplay", ::MyDisplay::Class_Version(), "MyDisplay.h", 48,
                  typeid(::MyDisplay), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MyDisplay::Dictionary, isa_proxy, 4,
                  sizeof(::MyDisplay) );
      instance.SetNew(&new_MyDisplay);
      instance.SetNewArray(&newArray_MyDisplay);
      instance.SetDelete(&delete_MyDisplay);
      instance.SetDeleteArray(&deleteArray_MyDisplay);
      instance.SetDestructor(&destruct_MyDisplay);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MyDisplay*)
   {
      return GenerateInitInstanceLocal(static_cast<::MyDisplay*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::MyDisplay*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr MyDisplay::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MyDisplay::Class_Name()
{
   return "MyDisplay";
}

//______________________________________________________________________________
const char *MyDisplay::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MyDisplay*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MyDisplay::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MyDisplay*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MyDisplay::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MyDisplay*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MyDisplay::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MyDisplay*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void MyDisplay::Streamer(TBuffer &R__b)
{
   // Stream an object of class MyDisplay.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MyDisplay::Class(),this);
   } else {
      R__b.WriteClassBuffer(MyDisplay::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MyDisplay(void *p) {
      return  p ? new(p) ::MyDisplay : new ::MyDisplay;
   }
   static void *newArray_MyDisplay(Long_t nElements, void *p) {
      return p ? new(p) ::MyDisplay[nElements] : new ::MyDisplay[nElements];
   }
   // Wrapper around operator delete
   static void delete_MyDisplay(void *p) {
      delete (static_cast<::MyDisplay*>(p));
   }
   static void deleteArray_MyDisplay(void *p) {
      delete [] (static_cast<::MyDisplay*>(p));
   }
   static void destruct_MyDisplay(void *p) {
      typedef ::MyDisplay current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::MyDisplay

namespace {
  void TriggerDictionaryInitialization_MyDisplayDict_Impl() {
    static const char* headers[] = {
"MyDisplay.h",
nullptr
    };
    static const char* includePaths[] = {
"/Users/ukose/sw/kits/root-install/include/",
"/Users/ukose/sw/Work/Magnet_a_la_babymind/MuonSpectrometerSim/Display/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MyDisplayDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$MyDisplay.h")))  MyDisplay;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MyDisplayDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "MyDisplay.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"MyDisplay", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MyDisplayDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MyDisplayDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MyDisplayDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MyDisplayDict() {
  TriggerDictionaryInitialization_MyDisplayDict_Impl();
}
