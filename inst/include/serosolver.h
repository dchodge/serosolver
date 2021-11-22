#include <RcppArmadillo.h>

using namespace Rcpp;

typedef void (func)(NumericVector &, const NumericVector &, const List &, const List &, const List &, const List &, const List &);
typedef void (*FunctionPtrAlt)(NumericVector &, const NumericVector &, const List &, const List &, const List &, const List &, const List &);
using abfunc = std::function<void(NumericVector &, const NumericVector &, const List &, const List &, const List &, const List &, const List &)>;
typedef void (* FunctionPtr)();

class FunctionPointer {
public:
    FunctionPtr ptr;
    FunctionPointer( FunctionPtr ptr_) : ptr(ptr_){}
} ;