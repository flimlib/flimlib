/* 
 * File:   LibraryWrapper.h
 * Author: aivar
 *
 * Created on August 27, 2010, 5:02 PM
 */

#ifndef _LIBRARYWRAPPER_H
#define	_LIBRARYWRAPPER_H

#ifdef	__cplusplus
extern "C" {
#endif

int markwardt_fit(double x_incr, double y[], int fit_start, int fit_end,
        double param[], int param_free[], int n_param);

#ifdef	__cplusplus
}
#endif

#endif	/* _LIBRARYWRAPPER_H */

