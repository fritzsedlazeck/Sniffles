// ***************************************************************************
// bamtools_variant.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011
// ---------------------------------------------------------------------------
// Provides a template-based variant type
// ---------------------------------------------------------------------------
// Modified from:
// variant_t - An Improved Variant Type Based on Member Templates
// (c) 2000 Fernando Cacciola
// Dr. Dobb's (http://www.ddj.com/cpp/184401293)
//
// * Modified to be in BamTools namespace, otherwise code is same. (DB)
// ***************************************************************************

#ifndef BAMTOOLS_VARIANT_H
#define BAMTOOLS_VARIANT_H

#include "utils/utils_global.h"
#include <stdexcept>
#include <string>
#include <typeinfo>

namespace BamTools {

class UTILS_EXPORT Variant {
  
    public:
        Variant(void) : data(NULL) { }
        
        Variant(const Variant& other) { 
            if ( other.data != NULL ) 
                other.data->AddRef();
            data = other.data;
        }

        ~Variant(void) { 
            if ( data != NULL ) 
                data->Release();
        }

        // NOTE: This code takes care of self-assignment.
        // DO NOT CHANGE THE ORDER of the statements.
        Variant& operator= (const Variant& rhs) {
            if ( rhs.data != NULL ) 
                rhs.data->AddRef();
            if ( data != NULL ) 
                data->Release();
            data = rhs.data;
            return *this;
        }

        // This member template constructor allows you to
        // instance a variant_t object with a value of any type.
        template<typename T>
        Variant(T v) 
            : data(new Impl<T>(v)) 
        { 
            data->AddRef(); 
        }

        // This generic conversion operator let you retrieve
        // the value held. To avoid template specialization conflicts,
        // it returns an instance of type T, which will be a COPY
        // of the value contained.
        template<typename T> 
        operator T() const { 
            return CastFromBase<T>(data)->data;
        }

        // This forms returns a REFERENCE and not a COPY, which
        // will be significant in some cases.
        template<typename T> 
        const T& get(void) const { 
            return CastFromBase<T>(data)->data; 
        }

        template<typename T> 
        bool is_type(void) const { 
            return typeid(*data)==typeid(Impl<T>); 
        }

        template<typename T> 
        bool is_type(T v) const { 
            return typeid(*data)==typeid(v); 
        }

    private:
        struct ImplBase {
                
            ImplBase() : refs(0) { }
            virtual ~ImplBase(void) { }
                
            void AddRef(void) { ++refs; }
            void Release(void) { 
                --refs;
                if ( refs == 0 ) delete this;
            }
                
            size_t refs;
        };

        template<typename T>
        struct Impl : ImplBase {
            Impl(T v) : data(v) { }
            ~Impl(void) { }
            T data;
        };

        // The following method is static because it doesn't
        // operate on variant_t instances.
        template<typename T> 
        static Impl<T>* CastFromBase(ImplBase* v) {
            // This upcast will fail if T is other than the T used
            // with the constructor of variant_t.
            Impl<T>* p = dynamic_cast< Impl<T>* > (v);
            if ( p == NULL ) 
                throw std::invalid_argument( typeid(T).name() + std::string(" is not a valid type") );
            return p;
        }

        ImplBase* data;
};

} // namespace BamTools

#endif // BAMTOOLS_VARIANT_H
