// -*- LSST-C++ -*-
#include "hpoly.h"

/******* Hermite methods ******************************************************/

double Hermite::returnHermiteN(int order, double x)
{
    //cout <<" order = "<< order <<endl;
    // coefficients are listed from zeroth order up to order n (n=order)
    vector<double> coeffs = getHermiteCoeffs(order);
    int ncoeffs = coeffs.size();
    if ( ncoeffs != (order+1) )
        throw ParmError("ERROR! Number of coefficient and order don't match");
    
    double Hn=0;
    for (int i=0; i<ncoeffs; i++) {
        double alpha = (double)i;
        Hn+=coeffs[i]*pow(x, alpha);
        }
    return Hn;

};


vector<double> Hermite::getHermiteCoeffs(int order)
{

    vector<double> coeffs;
    switch(order) {
        case 0: coeffs = c0_; return coeffs;
            break;
        case 1: coeffs = c1_; return coeffs;
            break;
        case 2: coeffs = c2_; return coeffs;
            break;  
        case 3: coeffs = c3_; return coeffs;
            break;
        case 4: coeffs = c4_; return coeffs;
            break;
        case 5: coeffs = c5_; return coeffs;
            break;
        case 6: coeffs = c6_; return coeffs;
            break;
        case 7: coeffs = c7_; return coeffs;
            break;
        case 8: coeffs = c8_; return coeffs;
            break;
        case 9: coeffs = c9_; return coeffs;
            break;
        case 10: coeffs = c10_; return coeffs;
            break;
        case 11: coeffs = c11_; return coeffs;
            break;
        case 12: coeffs = c12_; return coeffs;
            break;
        case 13: coeffs = c13_; return coeffs;
            break;
        case 14: coeffs = c14_; return coeffs;
            break;
        case 15: coeffs = c15_; return coeffs;
            break;
        case 16: coeffs = c16_; return coeffs;
            break;
        case 17: coeffs = c17_; return coeffs;
            break;
        case 18: coeffs = c18_; return coeffs;
            break; 
        case 19: coeffs = c19_; return coeffs;
            break;  
        case 20: coeffs = c20_; return coeffs;
            break;
        case 21: coeffs = c21_; return coeffs;
            break;
        case 22: coeffs = c22_; return coeffs;
            break;
        case 23: coeffs = c23_; return coeffs;
            break;
        case 24: coeffs = c24_; return coeffs;
            break;
        case 25: coeffs = c25_; return coeffs;
            break;
        default:
            coeffs = calculateHermiteCoeffs(order);
        }
};

vector<double> Hermite::calculateHermiteCoeffs(int order)
{
    // to do!
    cout <<"     Danger! shouldn't be in this function!"<<endl;
    if (order<=25)
        cout <<"     Warning! have pre-coded this order"<<endl;

    vector<double> coeffs(order+1);
    for (int i=0; i<=order/2; i++) {
        coeffs[order-2*i] = factorial(order)/(factorial(i)*factorial(order-2*i))*pow(2.,(order-2*i))*pow(-1,i);
        }
/*
endif else begin

  coeffs = lon64arr(n+1)
  coeffs[n] = 2LL^n ; This works up to (and including) n=30

  if n gt 24 then begin

    integers=ul64indgen(n+1)
    for i=1L,fix(n/2) do begin
      numerator=integers[(i>(n-2*i))+1:n]
      denominator=integers[1:i<(n-2*i)>1]
      for j=n_elements(denominator)-1,1,-1 do begin ; no need to check j=0 since that will always be 1
        multiple=where(numerator mod denominator[j] eq 0,n_multiple)
        if n_multiple gt 0 then begin
          numerator[(reverse(multiple))[0]]=numerator[(reverse(multiple))[0]]/denominator[j]
          denominator[j]=1LL
        endif
      endfor
      coeffs[n-2*i] = product(numerator,/INTEGER)/product(denominator,/INTEGER) * 2ULL^(n-2*i) * (-1)^i ; This works up to (and including) n=25. Above that, the final coefficient does not fit into a 64 bit signed integer (max is 2ULL^63).
    endfor

  endif else if n gt 20 then begin

    integers=ul64indgen(n+1)
    for i=1L,fix(n/2) do $ ; This works up to (and including) n=24, but the first term overflows above that
      coeffs[n-2*i] = product(integers[(i>(n-2*i))+1:n],/INTEGER) / factorial(i<(n-2*i)>1,/UL64) * 2ULL^(n-2*i) * (-1)^i
    ; Note that the ratio of the first two terms really does always give an integer. I don't have a rigorous mathematical proof of this, but it does.

  endif else begin

    for i=1L,fix(n/2) do $ ; This works up to (and including) n=20, but 21! doesn't fit into a UL64 integer
      coeffs[n-2*i] = factorial(n,/UL64) / ( factorial(i,/UL64) * factorial(n-2*i,/UL64) ) * 2ULL^(n-2*i) * (-1)^i

  endelse
endelse

return, coeffs*/

};

void Hermite::setLowOrderCoeffs()
{
    c0_.push_back(1);
    c1_.push_back(0); c1_.push_back(2);
    c2_.push_back(-2); c2_.push_back(0); c2_.push_back(4);
    c3_.push_back(0); c3_.push_back(-12); c3_.push_back(0); c3_.push_back(8); 
    c4_.push_back(12); c4_.push_back(0); c4_.push_back(-48); c4_.push_back(0); c4_.push_back(16);
    c5_.push_back(0); c5_.push_back(120); c5_.push_back(0); c5_.push_back(-160); c5_.push_back(0);
        c5_.push_back(32);
    c6_.push_back(-120); c6_.push_back(0); c6_.push_back(720); c6_.push_back(0); c6_.push_back(-480);
        c6_.push_back(0); c6_.push_back(64); 
    c7_.push_back(0); c7_.push_back(-1680); c7_.push_back(0); c7_.push_back(3360); c7_.push_back(0);
        c7_.push_back(-1344); c7_.push_back(0); c7_.push_back(128);
    c8_.push_back(1680); c8_.push_back(0); c8_.push_back(-13440); c8_.push_back(0); c8_.push_back(13440);
        c8_.push_back(0); c8_.push_back(-3584); c8_.push_back(0); c8_.push_back(256);
    c9_.push_back(0); c9_.push_back(30240); c9_.push_back(0); c9_.push_back(-80640); c9_.push_back(0); 
        c9_.push_back(48348); c9_.push_back(0); c9_.push_back(-9216); c9_.push_back(0); c9_.push_back(512);
    c10_.push_back(-30240); c10_.push_back(0); c10_.push_back(302400); c10_.push_back(0); c10_.push_back(-403200);
         c10_.push_back(0); c10_.push_back(161280); c10_.push_back(0); c10_.push_back(-23040); c10_.push_back(0);
            c10_.push_back(1024);
    c11_.push_back(0); c11_.push_back(-665280); c11_.push_back(0); c11_.push_back(2217600); c11_.push_back(0);
        c11_.push_back(-1774080); c11_.push_back(0); c11_.push_back(506880); c11_.push_back(0); c11_.push_back(-56320); 
            c11_.push_back(0); c11_.push_back(2048);
    c12_.push_back(665280); c12_.push_back(0); c12_.push_back(-7983360); c12_.push_back(0); c12_.push_back(13305600); 
        c12_.push_back(0); c12_.push_back(-7096320); c12_.push_back(0); c12_.push_back(1520640); c12_.push_back(0); 
            c12_.push_back(-135168); c12_.push_back(0); c12_.push_back(4096);
    c13_.push_back(0); c13_.push_back(17297280); c13_.push_back(0); c13_.push_back(-69189120); c13_.push_back(0); 
        c13_.push_back(69189120); c13_.push_back(0); c13_.push_back(-26357760); c13_.push_back(0); c13_.push_back(4392960); 
            c13_.push_back(0); c13_.push_back(-319488); c13_.push_back(0); c13_.push_back(8192);
    c14_.push_back(-17297280); c14_.push_back(0); c14_.push_back(242161920); c14_.push_back(0); c14_.push_back(-484323840); 
        c14_.push_back(0); c14_.push_back(322882560); c14_.push_back(0); c14_.push_back(-92252160); c14_.push_back(0); 
            c14_.push_back(12300288); c14_.push_back(0); c14_.push_back(-745472); c14_.push_back(0); c14_.push_back(16384);
    c15_.push_back(0); c15_.push_back(-518918400); c15_.push_back(0); c15_.push_back(2421619200); c15_.push_back(0); 
        c15_.push_back(-2905943040); c15_.push_back(0); c15_.push_back(1383782400); c15_.push_back(0); 
            c15_.push_back(-307507200); c15_.push_back(0); c15_.push_back(33546240); c15_.push_back(0); c15_.push_back(-1720320);
                c15_.push_back(0); c15_.push_back(32768);
    c16_.push_back(518918400); c16_.push_back(0); c16_.push_back(-8302694400); c16_.push_back(0); c16_.push_back(19372953600); 
        c16_.push_back(0); c16_.push_back(-15498362880); c16_.push_back(0); c16_.push_back(5535129600); c16_.push_back(0);
            c16_.push_back(-984023040); c16_.push_back(0); c16_.push_back(89456640); c16_.push_back(0); c16_.push_back(-3932160);
                c16_.push_back(0); c16_.push_back(65536);
    c17_.push_back(0); c17_.push_back(17643225600); c17_.push_back(0); c17_.push_back(-94097203200); c17_.push_back(0);
        c17_.push_back(131736084480); c17_.push_back(0); c17_.push_back(-75277762560); c17_.push_back(0); 
            c17_.push_back(20910489600); c17_.push_back(0); c17_.push_back(-3041525760); c17_.push_back(0); 
                c17_.push_back(233963520); c17_.push_back(0); c17_.push_back(-8912896); c17_.push_back(0); c17_.push_back(131072);
    c18_.push_back(-17643225600); c18_.push_back(0); c18_.push_back(317578060800); c18_.push_back(0); 
        c18_.push_back(-846874828800); c18_.push_back(0); c18_.push_back(790416506880); c18_.push_back(0); 
            c18_.push_back(-338749931520); c18_.push_back(0); c18_.push_back(75277762560); c18_.push_back(0); 
                c18_.push_back(-9124577280); c18_.push_back(0); c18_.push_back(601620480); c18_.push_back(0); 
                    c18_.push_back(-20054016); c18_.push_back(0); c18_.push_back(262144);
    c19_.push_back(0); c19_.push_back(-670442572800); c19_.push_back(0); c19_.push_back(4022655436800); c19_.push_back(0);
        c19_.push_back(-6436248698880); c19_.push_back(0); c19_.push_back(4290832465920); c19_.push_back(0); 
            c19_.push_back(-1430277488640); c19_.push_back(0); c19_.push_back(260050452480); c19_.push_back(0); 
                c19_.push_back(-26671841280); c19_.push_back(0); c19_.push_back(1524105216); c19_.push_back(0); 
                    c19_.push_back(-44826624); c19_.push_back(0); c19_.push_back(524288);
    c20_.push_back(670442572800); c20_.push_back(0); c20_.push_back(-13408851456000); c20_.push_back(0); 
        c20_.push_back(40226554368000); c20_.push_back(0); c20_.push_back(-42908324659200); c20_.push_back(0); 
            c20_.push_back(21454162329600); c20_.push_back(0); c20_.push_back(-5721109954560); c20_.push_back(0); 
                c20_.push_back(866834841600); c20_.push_back(0); c20_.push_back(-76205260800); c20_.push_back(0); 
                    c20_.push_back(3810263040); c20_.push_back(0); c20_.push_back(-99614720); c20_.push_back(0); 
                        c20_.push_back(1048576);
    c21_.push_back(0); c21_.push_back(28158588057600); c21_.push_back(0); c21_.push_back(-187723920384000); c21_.push_back(0);
        c21_.push_back(337903056691200); c21_.push_back(0); c21_.push_back(-257449947955200); c21_.push_back(0); 
            c21_.push_back(100119424204800); c21_.push_back(0); c21_.push_back(-21844238008320); c21_.push_back(0);
                c21_.push_back(2800543334400); c21_.push_back(0); c21_.push_back(-213374730240); c21_.push_back(0); 
                    c21_.push_back(9413591040); c21_.push_back(0); c21_.push_back(-220200960); c21_.push_back(0); 
                        c21_.push_back(2097152);
    c22_.push_back(-28158588057600); c22_.push_back(0); c22_.push_back(619488937267200); c22_.push_back(0); 
        c22_.push_back(-2064963124224000); c22_.push_back(0); c22_.push_back(2477955749068800); c22_.push_back(0);
            c22_.push_back(-1415974713753600); c22_.push_back(0); c22_.push_back(440525466501120); c22_.push_back(0);
                c22_.push_back(-80095539363840); c22_.push_back(0); c22_.push_back(8801707622400); c22_.push_back(0);
                    c22_.push_back(-586780508160); c22_.push_back(0); c22_.push_back(23011000320); c22_.push_back(0);       
                        c22_.push_back(-484442112); c22_.push_back(0); c22_.push_back(4194304);
    c23_.push_back(0); c23_.push_back(-1295295050649600); c23_.push_back(0); c23_.push_back(9498830371430400); 
        c23_.push_back(0); c23_.push_back(-18997660742860800); c23_.push_back(0); c23_.push_back(16283709208166400);
            c23_.push_back(0); c23_.push_back(-7237204092518400); c23_.push_back(0); c23_.push_back(1842197405368320); 
                c23_.push_back(0); c23_.push_back(-283414985441280); c23_.push_back(0); c23_.push_back(26991903375360); 
                    c23_.push_back(0); c23_.push_back(-1587759022080); c23_.push_back(0); c23_.push_back(55710842880); 
                        c23_.push_back(0); c23_.push_back(-1061158912); c23_.push_back(0); c23_.push_back(8388608);
    c24_.push_back(1295295050649600); c24_.push_back(0); c24_.push_back(-31087081215590400); c24_.push_back(0);
        c24_.push_back(113985964457164800); c24_.push_back(0); c24_.push_back(-151981285942886400); c24_.push_back(0);
            c24_.push_back(97702255248998400); c24_.push_back(0); c24_.push_back(-34738579644088320); c24_.push_back(0);
                c24_.push_back(7368789621473280); c24_.push_back(0); c24_.push_back(-971708521512960); c24_.push_back(0);
                    c24_.push_back(80975710126080); c24_.push_back(0); c24_.push_back(-4234024058880); c24_.push_back(0);
                        c24_.push_back(133706022912); c24_.push_back(0); c24_.push_back(-2315255808); c24_.push_back(0); 
                            c24_.push_back(16777216);
    c25_.push_back(0); c25_.push_back(64764752532480000); c25_.push_back(0); c25_.push_back(-518118020259840000); 
        c25_.push_back(0); c25_.push_back(1139859644571648000); c25_.push_back(0); c25_.push_back(-1085580613877760000); 
            c25_.push_back(0); c25_.push_back(542790306938880000); c25_.push_back(0); c25_.push_back(-157902634745856000); 
                c25_.push_back(0); c25_.push_back(28341498544128000); c25_.push_back(0); c25_.push_back(-3239028405043200); 
                    c25_.push_back(0); c25_.push_back(238163853312000); c25_.push_back(0); c25_.push_back(-11142168576000); 
                        c25_.push_back(0); c25_.push_back(318347673600); c25_.push_back(0); c25_.push_back(-5033164800); 
                            c25_.push_back(0); c25_.push_back(33554432);

};
