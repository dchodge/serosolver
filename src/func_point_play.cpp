#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]
typedef void (* ab_kin_type)(int x);
//typedef vec (*funcPtr)(const vec& x);

// [[Rcpp::export]]
void fun_cpp(int x) {
    x = x + 3;
}
/*
// [[Rcpp::export]]
void callViaXPtr(const int x, SEXP xpsexp) {
    XPtr<ab_kin_type> xpfun(xpsexp);
    ab_kin_type fun = *xpfun;
    fun(x);
}

// [[Rcpp::export]]
XPtr<ab_kin_type> putFunPtrInXPtr(ab_kin_type func) {
    return(XPtr<ab_kin_type>(new funcPtr(&func)));
}
*/
 /*
class AbFuncs
{
public:
    AbFuncs(double x) {}

    std::function<void(int)> fcn;

    double random = 5;
};

RCPP_MODULE(AbFuncsMod) {
    class_<AbFuncs>( "AbFuncs" )

    .constructor<double >()
    .field("random", &AbFuncs::random)
    .method( "fcn", &AbFuncs::fcn)
    ;
}
*/


// define the structure
struct xptr_data{
  SEXP xptr;   // SEXP S-expression, a common datatype for all R objects
};

// a minimal function (user-defined)
// [[Rcpp::export]]
void timesTwo(NumericVector x) {
   x = x * 2;
}

// pointer to function defined
//typedef NumericVector (*funcPtr) (NumericVector y);

// [[Rcpp::export]]
SEXP makeFuncPtr(SEXP xpsexp) {
  auto fcnPtr { &xpsexp }; 
  return *fcnPtr;
}

///void evalPtr(const int x, SEXP funcPtr) {
   // funcPtr(x);
//}

/*
// this function will be in the package
NumericVector call_by_xptr_struct(NumericVector y, void* user_data){

  struct xptr_data *my_rhs_ptr = (struct xptr_data*)user_data;
  SEXP xpsexp = (*my_rhs_ptr).xptr;

  // use function pointer to get the derivatives
  XPtr<funcPtr> rhs_xptr(xpsexp);
  funcPtr rhs_fun = *rhs_xptr;

  // use the function to calculate value of RHS ----
  return(rhs_fun(y));
}


// using xptr to evaluate function - this will be exported
// from the package
//[[Rcpp::export]]
NumericVector xptr_call_struct(NumericVector y, SEXP xpsexp){

  struct xptr_data my_xptr = {NULL};

  my_xptr.xptr = xpsexp;
  return call_by_xptr_struct(y, (void*)&my_xptr);
}

*/