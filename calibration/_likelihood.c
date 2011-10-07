#include <Python.h>
#include <numpy/arrayobject.h>

#ifdef USEOPENMP
#include <omp.h>
#endif


/* ====================== //
// PYTHON INITIALIZATIONS //
// ====================== */

PyMODINIT_FUNC init_likelihood(void);
static PyObject *likelihood_lnlikelihood(PyObject *self, PyObject *args);
static PyObject *likelihood_lnoddsvar(PyObject *self, PyObject *args);
static PyObject *likelihood_lnoddsbad(PyObject *self, PyObject *args);
static PyObject *likelihood_lnlikeratiobad(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"lnlikelihood", likelihood_lnlikelihood, METH_VARARGS,"Calculate the lnlikelihood."},
    {"lnoddsvar", likelihood_lnoddsvar, METH_VARARGS,"Calculate the odds that a star is variable."},
    {"lnoddsbad", likelihood_lnoddsbad, METH_VARARGS,"Calculate the odds that a specific measurement is bad."},
    {"lnlikeratiobad", *likelihood_lnlikeratiobad, METH_VARARGS,"Calculate the likelihood ratio that a specific measurement is bad."},
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

// ================ //
// PYTHON INTERFACE //
// ================ //

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

/* ================ //
// TYPE DEFINITIONS //
// ================ */

typedef struct patchmodel {
    PyObject *df, *di;
    double *flux, *ivar;
    int nstars, nobs;
    PyObject *mz, *mf;
    double *zero, *fstar;
    double jitterabs2, jitterrel2, Q2, Pvar, sigbad2, Pbad;
} PatchModel;

PatchModel *PatchModel_init(PyObject *model)
{
    PatchModel *self = malloc(sizeof(PatchModel));

    // parse the data/model classes
    self->df = _fromNPYArray(model,"_f");
    self->di = _fromNPYArray(model,"_ivar_f");
    if (self->df == NULL || self->di == NULL)  {
        return NULL;
    }

    self->flux = (double *)PyArray_DATA(self->df);
    self->ivar = (double *)PyArray_DATA(self->di);
    self->nobs = PyArray_DIMS(self->df)[0];
    self->nstars = PyArray_DIMS(self->df)[1];

    // parse the data/model classes
    self->mz = _fromNPYArray(model,"_f0");
    self->mf = _fromNPYArray(model,"_fstar");
    if (self->mz == NULL || self->mf == NULL)
        return NULL;

    self->zero = (double *)PyArray_DATA(self->mz);
    self->fstar = (double *)PyArray_DATA(self->mf);

    // model parameters
    self->jitterabs2 = _fromPyDouble(model,"_jitterabs2");
    self->jitterrel2 = _fromPyDouble(model,"_jitterrel2");
    self->Q2         = _fromPyDouble(model,"_Q2");
    self->Pvar       = _fromPyDouble(model,"_pvar");
    self->sigbad2    = _fromPyDouble(model,"_sigbad2");
    self->Pbad       = _fromPyDouble(model,"_pbad");

    return self;
}

void PatchModel_destroy(PatchModel *self)
{
    if (self != NULL) {
        Py_XDECREF(self->df);
        Py_XDECREF(self->di);

        Py_XDECREF(self->mz);
        Py_XDECREF(self->mf);

        free(self);
    }
}

// ============== //
// MATH FUNCTIONS //
// ============== //

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

// ======================= //
// LIKELIHOOD CALCULATIONS //
// ======================= //

void get_lnpgood_and_lnpbad_and_lnpvargood(int i, int alpha,
        double *lnpgood, double *lnpbad, double *lnpvargood,
        PatchModel *model)
{
    int ind = i*model->nstars+alpha; // data is TRANSPOSED from what we want!
    *lnpgood = 0.0; *lnpbad = 0.0; *lnpvargood = 0.0;

    if (model->ivar[ind] > 0.0) {
        double ff = model->zero[i]*model->fstar[alpha]; // outer(zero,flux)
        double sig2   = 1.0/model->ivar[ind];
        double delta2 = model->jitterabs2 + model->jitterrel2*ff*ff;

        *lnpgood    = _lnnormal(model->flux[ind],ff,sig2+delta2);
        *lnpbad     = _lnnormal(model->flux[ind],ff,sig2+delta2+model->sigbad2);
        *lnpvargood = _lnnormal(model->flux[ind],ff,sig2+delta2+model->Q2*ff*ff);
    }
}

void get_lnpvar_and_lnpconst(int alpha, double *lnpvar, double *lnpconst,
        PatchModel *model)
{
    int i;
    *lnpconst = 0.0;
    *lnpvar   = 0.0;
    for (i = 0; i < model->nobs; i++) {
        double lnpgood,lnpbad,lnpvargood;
        get_lnpgood_and_lnpbad_and_lnpvargood(i,alpha,
                &lnpgood,&lnpbad,&lnpvargood,
                model);
        *lnpconst += _logsumexp(log(1-model->Pbad)+lnpgood,
                                 log(model->Pbad)+lnpbad);
        *lnpvar   += _logsumexp(log(1-model->Pbad)+lnpvargood,
                                 log(model->Pbad)+lnpbad);
    }
}

double lnlikelihood(PatchModel *model)
{
    int alpha;
    double lnlike = 0.0;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(alpha) \
    schedule(static) reduction(+:lnlike)
#endif
    for (alpha = 0; alpha < model->nstars; alpha++) {
        double lnpconst, lnpvar;
        get_lnpvar_and_lnpconst(alpha, &lnpvar, &lnpconst, model);
        lnlike += _logsumexp(log(1-model->Pvar)+lnpconst,log(model->Pvar)+lnpvar);
    }
    return lnlike;
}

void get_lnoddsvar(PatchModel *model, double *lnoddsvar)
{
    int alpha;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(alpha) \
    schedule(static)
#endif
    for (alpha = 0; alpha < model->nstars; alpha++) {
        double lnpconst, lnpvar;

        get_lnpvar_and_lnpconst(alpha, &lnpvar, &lnpconst, model);

        lnoddsvar[alpha] = log(model->Pvar)-log(1-model->Pvar)
            + lnpvar - lnpconst;
    }
}

void get_lnlikeratiobad(PatchModel *model, double *lnlikeratiobad)
{
    int i,alpha;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(alpha) \
    schedule(static)
#endif
    for (i = 0; i < model->nobs; i++) {
        for (alpha = 0; alpha < model->nstars; alpha++) {
            int ind = i*model->nstars+alpha;

            double lnpgood,lnpbad,lnpvargood;
            double lnpconst,lnpvar,lntotgood;

            get_lnpgood_and_lnpbad_and_lnpvargood(i,alpha,
                    &lnpgood,&lnpbad,&lnpvargood, model);
            get_lnpvar_and_lnpconst(alpha, &lnpvar, &lnpconst, model);

            //lntotgood = _logsumexp(lnpvar+lnpvargood,lnpconst+lnpgood);
            lntotgood = lnpgood;
            lnlikeratiobad[ind] = lnpbad - lntotgood;
        }
    }
}

void get_lnoddsbad(PatchModel *model, double *lnoddsbad)
{
    int i,alpha;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(alpha) \
    schedule(static)
#endif
    for (i = 0; i < model->nobs; i++) {
        for (alpha = 0; alpha < model->nstars; alpha++) {
            int ind = i*model->nstars+alpha;

            double lnpgood,lnpbad,lnpvargood;
            double lnpconst,lnpvar,lntotgood, lntotbad;

            get_lnpgood_and_lnpbad_and_lnpvargood(i,alpha,
                    &lnpgood,&lnpbad,&lnpvargood, model);
            get_lnpvar_and_lnpconst(alpha, &lnpvar, &lnpconst, model);

            lntotgood = _logsumexp(log(model->Pvar)+lnpvar+lnpvargood,
                                lnpconst+lnpgood+log(1-model->Pvar));
            lntotbad = _logsumexp(log(model->Pvar)+lnpvar+lnpbad,
                                lnpconst+lnpbad+log(1-model->Pvar));

            lnoddsbad[ind] = log(model->Pbad)-log(1-model->Pbad)
                + lntotbad - lntotgood;
        }
    }
}

// ============== //
// MODULE METHODS //
// ============== //

static PyObject *likelihood_lnlikelihood(PyObject *self, PyObject *args)
{
    // 2011-07-06 - Created by Dan F-M
    PyObject *model0 = NULL;
    PatchModel *model = NULL;
    if (!PyArg_ParseTuple(args, "O", &model0))
        goto fail;

    model = PatchModel_init(model0);
    if (model == NULL)
        goto fail;

    // Calcuate the likelihood
    double lnlike = lnlikelihood(model);
    PyObject *result = PyFloat_FromDouble(lnlike);

    // clean up!
    PatchModel_destroy(model);
    return result;

fail:
    PyErr_SetString(PyExc_RuntimeError, "Incorrect input for likelihood");
    
    // clean up and fail
    PatchModel_destroy(model);
    
    return NULL;
}

static PyObject *likelihood_lnoddsvar(PyObject *self, PyObject *args)
{
    PyObject *model0 = NULL, *lnoddsvar_obj = NULL;
    PyObject *lnoddsvar = NULL;
    PatchModel *model = NULL;    
    if (!PyArg_ParseTuple(args, "OO", &model0, &lnoddsvar_obj))
        goto fail;
    
    model = PatchModel_init(model0);
    if (model == NULL)
        goto fail;

    lnoddsvar = PyArray_FROM_OTF(lnoddsvar_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    if (lnoddsvar == NULL)
        goto fail;
    get_lnoddsvar(model, (double *)PyArray_DATA(lnoddsvar));

    // clean up!
    Py_DECREF(lnoddsvar);
    PatchModel_destroy(model);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    PyErr_SetString(PyExc_RuntimeError, "Incorrect input for lnoddsvar");
    
    Py_XDECREF(lnoddsvar);
    PatchModel_destroy(model);
    
    return NULL;
}

static PyObject *likelihood_lnoddsbad(PyObject *self, PyObject *args)
{
    PyObject *model0 = NULL, *lnoddsbad_obj = NULL, *lnoddsbad = NULL;
    PatchModel *model = NULL;
    if (!PyArg_ParseTuple(args, "OO", &model0, &lnoddsbad_obj))
        goto fail;
    
    model = PatchModel_init(model0);
    if (model == NULL)
        goto fail;

    lnoddsbad = PyArray_FROM_OTF(lnoddsbad_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    if (lnoddsbad == NULL)
        goto fail;
    get_lnoddsbad(model, (double *)PyArray_DATA(lnoddsbad));

    // clean up!
    Py_DECREF(lnoddsbad);
    PatchModel_destroy(model);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    PyErr_SetString(PyExc_RuntimeError, "Incorrect input for lnoddsbad");
    
    Py_XDECREF(lnoddsbad);
    PatchModel_destroy(model);
    
    return NULL;
}

static PyObject *likelihood_lnlikeratiobad(PyObject *self, PyObject *args)
{
    PyObject *model0 = NULL, *lnlikeratiobad_obj = NULL, *lnlikeratiobad = NULL;
    PatchModel *model = NULL;
    if (!PyArg_ParseTuple(args, "OO", &model0, &lnlikeratiobad_obj))
        goto fail;
    
    model = PatchModel_init(model0);
    if (model == NULL)
        goto fail;

    lnlikeratiobad = PyArray_FROM_OTF(lnlikeratiobad_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    if (lnlikeratiobad == NULL)
        goto fail;
    get_lnlikeratiobad(model, (double *)PyArray_DATA(lnlikeratiobad));
    
    // clean up!
    Py_DECREF(lnlikeratiobad);
    PatchModel_destroy(model);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    PyErr_SetString(PyExc_RuntimeError, "Incorrect input for lnlikeratiobad");
    
    Py_XDECREF(lnlikeratiobad);
    PatchModel_destroy(model);
    
    return NULL;
}

