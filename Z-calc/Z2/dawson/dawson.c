/* specfunc/dawson.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/**
 * ИЗМЕНЕНИЯ: 
 * Удалено все, что связано с кодами ошибок.
 * Добавлен внутрь cheb_eval.c.
 * Добавлен #include <math.h> и M_SQRT2.
 */

#include "gsl_machine.h"
#include "gsl_sf_dawson.h"

#include "chebyshev.h"
#include "error.h"

#include <math.h>

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872421 /* sqrt(2) */
#endif

/* Based on ddaws.f, Fullerton, W., (LANL) */

/* Из cheb_eval.c. */
inline static void
cheb_eval_e(const cheb_series *cs,
            const double x,
            gsl_sf_result *result) {
    int j;
    double d  = 0.0;
    double dd = 0.0;

    double y  = (2.0 * x - cs->a - cs->b) / (cs->b - cs->a);
    double y2 = 2.0 * y;

    double e = 0.0;

    for (j = cs->order; j >= 1; j--) {
        double temp = d;
        d           = y2 * d - dd + cs->c[j];
        e += fabs(y2 * temp) + fabs(dd) + fabs(cs->c[j]);
        dd = temp;
    }

    {
        double temp = d;
        d           = y * d - dd + 0.5 * cs->c[0];
        e += fabs(y * temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
    }

    result->val = d;
    result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);
}

/* Chebyshev expansions

 Series for DAW        on the interval  0.          to  1.00000E+00
                                        with weighted error   8.95E-32
                                         log weighted error  31.05
                               significant figures required  30.41
                                    decimal places required  31.71

 Series for DAW2       on the interval  0.          to  1.60000E+01
                                        with weighted error   1.61E-32
                                         log weighted error  31.79
                               significant figures required  31.40
                                    decimal places required  32.62

 Series for DAWA       on the interval  0.          to  6.25000E-02
                                        with weighted error   1.97E-32
                                         log weighted error  31.71
                               significant figures required  29.79
                                    decimal places required  32.64
*/
static double daw_data[21] = {
    -0.6351734375145949201065127736293e-02,
    -0.2294071479677386939899824125866e+00,
    0.2213050093908476441683979161786e-01,
    -0.1549265453892985046743057753375e-02,
    0.8497327715684917456777542948066e-04,
    -0.3828266270972014924994099521309e-05,
    0.1462854806250163197757148949539e-06,
    -0.4851982381825991798846715425114e-08,
    0.1421463577759139790347568183304e-09,
    -0.3728836087920596525335493054088e-11,
    0.8854942961778203370194565231369e-13,
    -0.1920757131350206355421648417493e-14,
    0.3834325867246327588241074439253e-16,
    -0.7089154168175881633584099327999e-18,
    0.1220552135889457674416901120000e-19,
    -0.1966204826605348760299451733333e-21,
    0.2975845541376597189113173333333e-23,
    -0.4247069514800596951039999999999e-25,
    0.5734270767391742798506666666666e-27,
    -0.7345836823178450261333333333333e-29,
    0.8951937667516552533333333333333e-31
};
static cheb_series daw_cs = {
    daw_data,
    15, /* 20, */
    -1, 1,
    9
};

static double daw2_data[45] = {
    -0.56886544105215527114160533733674e-01,
    -0.31811346996168131279322878048822e+00,
    0.20873845413642236789741580198858e+00,
    -0.12475409913779131214073498314784e+00,
    0.67869305186676777092847516423676e-01,
    -0.33659144895270939503068230966587e-01,
    0.15260781271987971743682460381640e-01,
    -0.63483709625962148230586094788535e-02,
    0.24326740920748520596865966109343e-02,
    -0.86219541491065032038526983549637e-03,
    0.28376573336321625302857636538295e-03,
    -0.87057549874170423699396581464335e-04,
    0.24986849985481658331800044137276e-04,
    -0.67319286764160294344603050339520e-05,
    0.17078578785573543710504524047844e-05,
    -0.40917551226475381271896592490038e-06,
    0.92828292216755773260751785312273e-07,
    -0.19991403610147617829845096332198e-07,
    0.40963490644082195241210487868917e-08,
    -0.80032409540993168075706781753561e-09,
    0.14938503128761465059143225550110e-09,
    -0.26687999885622329284924651063339e-10,
    0.45712216985159458151405617724103e-11,
    -0.75187305222043565872243727326771e-12,
    0.11893100052629681879029828987302e-12,
    -0.18116907933852346973490318263084e-13,
    0.26611733684358969193001612199626e-14,
    -0.37738863052129419795444109905930e-15,
    0.51727953789087172679680082229329e-16,
    -0.68603684084077500979419564670102e-17,
    0.88123751354161071806469337321745e-18,
    -0.10974248249996606292106299624652e-18,
    0.13261199326367178513595545891635e-19,
    -0.15562732768137380785488776571562e-20,
    0.17751425583655720607833415570773e-21,
    -0.19695006967006578384953608765439e-22,
    0.21270074896998699661924010120533e-23,
    -0.22375398124627973794182113962666e-24,
    0.22942768578582348946971383125333e-25,
    -0.22943788846552928693329592319999e-26,
    0.22391702100592453618342297600000e-27,
    -0.21338230616608897703678225066666e-28,
    0.19866196585123531518028458666666e-29,
    -0.18079295866694391771955199999999e-30,
    0.16090686015283030305450666666666e-31
};
static cheb_series daw2_cs = {
    daw2_data,
    32, /* 44, */
    -1, 1,
    21
};

static double dawa_data[75] = {
    0.1690485637765703755422637438849e-01,
    0.8683252278406957990536107850768e-02,
    0.2424864042417715453277703459889e-03,
    0.1261182399572690001651949240377e-04,
    0.1066453314636176955705691125906e-05,
    0.1358159794790727611348424505728e-06,
    0.2171042356577298398904312744743e-07,
    0.2867010501805295270343676804813e-08,
    -0.1901336393035820112282492378024e-09,
    -0.3097780484395201125532065774268e-09,
    -0.1029414876057509247398132286413e-09,
    -0.6260356459459576150417587283121e-11,
    0.8563132497446451216262303166276e-11,
    0.3033045148075659292976266276257e-11,
    -0.2523618306809291372630886938826e-12,
    -0.4210604795440664513175461934510e-12,
    -0.4431140826646238312143429452036e-13,
    0.4911210272841205205940037065117e-13,
    0.1235856242283903407076477954739e-13,
    -0.5788733199016569246955765071069e-14,
    -0.2282723294807358620978183957030e-14,
    0.7637149411014126476312362917590e-15,
    0.3851546883566811728777594002095e-15,
    -0.1199932056928290592803237283045e-15,
    -0.6313439150094572347334270285250e-16,
    0.2239559965972975375254912790237e-16,
    0.9987925830076495995132891200749e-17,
    -0.4681068274322495334536246507252e-17,
    -0.1436303644349721337241628751534e-17,
    0.1020822731410541112977908032130e-17,
    0.1538908873136092072837389822372e-18,
    -0.2189157877645793888894790926056e-18,
    0.2156879197938651750392359152517e-20,
    0.4370219827442449851134792557395e-19,
    -0.8234581460977207241098927905177e-20,
    -0.7498648721256466222903202835420e-20,
    0.3282536720735671610957612930039e-20,
    0.8858064309503921116076561515151e-21,
    -0.9185087111727002988094460531485e-21,
    0.2978962223788748988314166045791e-22,
    0.1972132136618471883159505468041e-21,
    -0.5974775596362906638089584995117e-22,
    -0.2834410031503850965443825182441e-22,
    0.2209560791131554514777150489012e-22,
    -0.5439955741897144300079480307711e-25,
    -0.5213549243294848668017136696470e-23,
    0.1702350556813114199065671499076e-23,
    0.6917400860836148343022185660197e-24,
    -0.6540941793002752512239445125802e-24,
    0.6093576580439328960371824654636e-25,
    0.1408070432905187461501945080272e-24,
    -0.6785886121054846331167674943755e-25,
    -0.9799732036214295711741583102225e-26,
    0.2121244903099041332598960939160e-25,
    -0.5954455022548790938238802154487e-26,
    -0.3093088861875470177838847232049e-26,
    0.2854389216344524682400691986104e-26,
    -0.3951289447379305566023477271811e-27,
    -0.5906000648607628478116840894453e-27,
    0.3670236964668687003647889980609e-27,
    -0.4839958238042276256598303038941e-29,
    -0.9799265984210443869597404017022e-28,
    0.4684773732612130606158908804300e-28,
    0.5030877696993461051647667603155e-29,
    -0.1547395051706028239247552068295e-28,
    0.6112180185086419243976005662714e-29,
    0.1357913399124811650343602736158e-29,
    -0.2417687752768673088385304299044e-29,
    0.8369074582074298945292887587291e-30,
    0.2665413042788979165838319401566e-30,
    -0.3811653692354890336935691003712e-30,
    0.1230054721884951464371706872585e-30,
    0.4622506399041493508805536929983e-31,
    -0.6120087296881677722911435593001e-31,
    0.1966024640193164686956230217896e-31
};
static cheb_series dawa_cs = {
    dawa_data,
    34, /* 74, */
    -1, 1,
    12
};

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

void gsl_sf_dawson_e(double x, gsl_sf_result *result) {
    const double xsml = 1.225 * GSL_SQRT_DBL_EPSILON;
    const double xbig = 1.0 / (M_SQRT2 * GSL_SQRT_DBL_EPSILON);
    const double xmax = 0.1 * GSL_DBL_MAX;

    const double y = fabs(x);

    if (y < xsml) {
        result->val = x;
        result->err = 0.0;
        // return GSL_SUCCESS;
    } else if (y < 1.0) {
        gsl_sf_result result_c;
        cheb_eval_e(&daw_cs, 2.0 * y * y - 1.0, &result_c);
        result->val = x * (0.75 + result_c.val);
        result->err = y * result_c.err;
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        // return GSL_SUCCESS;
    } else if (y < 4.0) {
        gsl_sf_result result_c;
        cheb_eval_e(&daw2_cs, 0.125 * y * y - 1.0, &result_c);
        result->val = x * (0.25 + result_c.val);
        result->err = y * result_c.err;
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        // return GSL_SUCCESS;
    } else if (y < xbig) {
        gsl_sf_result result_c;
        cheb_eval_e(&dawa_cs, 32.0 / (y * y) - 1.0, &result_c);
        result->val = (0.5 + result_c.val) / x;
        result->err = result_c.err / y;
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        // return GSL_SUCCESS;
    } else if (y < xmax) {
        result->val = 0.5 / x;
        result->err = 2.0 * GSL_DBL_EPSILON * result->val;
        // return GSL_SUCCESS;
    }
    // else {
    //   UNDERFLOW_ERROR(result);
    // }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_dawson(double x) {
    gsl_sf_result result;
    gsl_sf_dawson_e(x, &result);
    return result.val;
}
