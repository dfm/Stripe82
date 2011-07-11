#include <Python.h>
#include <numpy/arrayobject.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

// ====================== //
// PYTHON INITIALIZATIONS //
// ====================== //

PyMODINIT_FUNC init_likelihood(void);
static PyObject *likelihood_lnlikelihood(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"lnlikelihood", likelihood_lnlikelihood, METH_VARARGS,"Calculate the lnlikelihood."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_likelihood(void)
{
    PyObject *m;

    m = Py_InitModule("_likelihood", module_methods);
    if (m == NULL)
        return;

    import_array();
}

static const double lisqrt2pi = -0.91893853320467267; // -0.5*log(2.0*M_PI);
double _lnnormal(double x, double mu, double var)
{
    double diff = x-mu;
    return -0.5*diff*diff/var - 0.5*log(var) + lisqrt2pi;
}
double _logsumexp(double a, double b)
{
    if (a > b)
        return a + log(1+exp(b-a));
    return b + log(1+exp(a-b));
}
double _fromPyDouble(PyObject *obj, char *attr)
{
    PyObject *number = PyObject_GetAttrString(obj,attr);
    double result = PyFloat_AsDouble(number);
    Py_DECREF(number);
    return result;
}
PyObject *_fromNPYArray(PyObject *obj, char *attr)
{
    PyObject *arr = PyObject_GetAttrString(obj,attr);
    PyObject *result = PyArray_FROM_OTF(arr, NPY_DOUBLE, NPY_IN_ARRAY);
    Py_DECREF(arr);
    return result;
}

static PyObject *likelihood_lnlikelihood(PyObject *self, PyObject *args)
{
    // 2011-07-06 - Created by Dan F-M
    //
    // WARNING: For speed, the data and model arrays are actually transposed
    //          with respect to each other ... be cautious!

    /* ================== *
     * Python MUMBO JUMBO *
     * ================== */

    PyObject *model = NULL, *data = NULL;
    if (!PyArg_ParseTuple(args, "OO", &model, &data)) {
        PyErr_SetString(PyExc_RuntimeError,
                "Incorrect input for likelihood");
        return NULL;
    }

    // parse the data/model classes
    PyObject *dataflux = _fromNPYArray(data,"flux");
    PyObject *dataivar = _fromNPYArray(data,"ivar");
    PyObject *mz = _fromNPYArray(model,"zero");
    PyObject *mf = _fromNPYArray(model,"flux");
    if (dataflux == NULL || dataivar == NULL || mz == NULL || mf == NULL)
        goto fail;

    // model parameters
    double jitterabs2 = _fromPyDouble(model,"jitterabs2");
    double jitterrel2 = _fromPyDouble(model,"jitterrel2");
    double sigvar2    = _fromPyDouble(model,"sigvar2");
    double Pvar       = _fromPyDouble(model,"pvar");
    double sigbad2    = _fromPyDouble(model,"sigbad2");
    double Pbad       = _fromPyDouble(model,"pbad");
    
    /* =================================== *
     * Do the calculations in C for speed! *
     * =================================== */

    int i,j;    
    double *mzero = (double *)PyArray_DATA(mz);
    double *mflux = (double *)PyArray_DATA(mf);
    double *df = (double *)PyArray_DATA(dataflux);
    double *ivar = (double *)PyArray_DATA(dataivar);
    int nobs = PyArray_DIMS(mz)[0], nstars = PyArray_DIMS(mf)[0];

    // calculate lnprob1
    double lnprob = 0.0;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(i,j) \
    schedule(static) reduction(+:lnprob)
#endif
    for (i = 0; i < nstars; i++) {
        double lnpconst = 0.0, lnpvar = 0.0;
        for (j = 0; j < nobs; j++) {
            int ind = j*nstars+i; // data is TRANSPOSED from what we want!
            if (ivar[ind] > 0.0) {
                double ff = mzero[j]*mflux[i]; // outer(zero,flux)
                double sig2   = 1.0/ivar[ind];
                double delta2 = jitterabs2 + jitterrel2*ff*ff;

                double lnpgood    = _lnnormal(df[ind],ff,sig2+delta2);
                double lnpbad     = _lnnormal(df[ind],ff,sig2+delta2+sigbad2);
                double lnpvargood = _lnnormal(df[ind],ff,sig2+delta2+sigvar2);

                lnpconst += _logsumexp(log(1-Pbad)+lnpgood,
                                         log(Pbad)+lnpbad);
                lnpvar   += _logsumexp(log(1-Pbad)+lnpvargood,
                                         log(Pbad)+lnpbad);
            }
        }
        lnprob = lnprob + _logsumexp(log(1-Pvar)+lnpconst,log(Pvar)+lnpvar);
    }
    PyObject *result = PyFloat_FromDouble(lnprob);
    //Py_INCREF(result);

    // clean up!
    Py_DECREF(dataflux);
    Py_DECREF(dataivar);
    Py_DECREF(mz);
    Py_DECREF(mf);

    return result;

fail:
    // clean up and fail
    Py_XDECREF(dataflux);
    Py_XDECREF(dataivar);
    Py_XDECREF(mz);
    Py_XDECREF(mf);
    
    PyErr_SetString(PyExc_RuntimeError,
                "Likelihood calculation failed");
    return NULL;
}




