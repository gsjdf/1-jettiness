#ifndef __lha_pdf_h__
#define __lha_pdf_h__ 1


#include <bits/dis-process.h>
#include "LHAPDF/LHAPDF.h"


//----- used namespaces -----
using namespace nlo;
using namespace std;


class lha_pdf
  : public pdf_and_coupling_dis
{
public:
  //   constructor
  explicit lha_pdf(string name,LHAPDF::SetType type,unsigned int mem = 0)
    { //LHAPDF::initPDFSet("cteq61", LHAPDF::LHPDF, mem);
      LHAPDF::initPDFSet(name, type, mem);}

    explicit lha_pdf()
    { LHAPDF::initPDFSet("NNPDF30_nlo_as_0118", LHAPDF::LHGRID, 0);} 
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
// test of fixed alphas
    return LHAPDF::alphasPDF(std::sqrt(mr2))/6.28318530717958647692;
  }
  
  //   the parton distribution function
  void hadron(double x, double Q2, unsigned int, unsigned int, double *f) {
    vector<double> __f = LHAPDF::xfx(x, sqrt(Q2));
    for(int i=-6; i <= 6; i++) f[i] = __f[6+i]/x;

    //f[0] = 0.;
    //f[-5] = 0.;
    //f[-4] = 0.;
    //f[-3] = 0.;
    //f[-2] = 0.;
    //f[-1] = 0.;
    //f[5] = 0.;
    //f[4] = 0.;
    //f[3] = 0.;
    ////f[2] = 0.;
    //f[1] = 0.;
  }
  
private:

};



#endif
