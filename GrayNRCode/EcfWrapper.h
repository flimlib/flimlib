/* 
 * File:   EcfWrapper.h
 * Author: Aivar Grislis
 *
 * Created on September 3, 2010, 5:21 PM
 */

#ifndef _ECFWRAPPER_H
#define	_ECFWRAPPER_H

#ifdef	__cplusplus
extern "C" {
#endif

int RLD_fit(
        double x_inc,
        double y[],
        int fit_start,
        int fit_end,
        double instr[],
        int n_instr,
        double sig[],
        double *z,
        double *a,
        double *tau,
        double fitted[],
        double *chi_square,
        double chi_square_target
        );

int LMA_fit(
        double x_inc,
        double y[],
        int fit_start,
        int fit_end,
        double instr[],
        int n_instr,
        double sig[],
        double param[],
        int param_free[],
        int n_param,
        double fitted[],
        double *chi_square,
        double chi_square_target
        );

#ifdef	__cplusplus
}
#endif

#endif	/* _ECFWRAPPER_H */

