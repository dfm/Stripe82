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

typedef struct photodata {
    PyObject *df, *di;
    double *flux, *ivar;
    int nstars, nobs;
} PhotoData;

typedef struct photomodel {
    PyObject *mz, *mf;
    double *zero, *flux;
    double jitterabs2, jitterrel2, Q2, Pvar, sigbad2, Pbad;
} PhotoModel;

PhotoData *PhotoData_init(PyObject *data)
{
    PhotoData *self = malloc(sizeof(PhotoData));

    // parse the data/model classes
    self->df = _fromNPYArray(data,"flux");
    self->di = _fromNPYArray(data,"ivar");
    if (self->df == NULL || self->di == NULL)  {
        return NULL;
    }

    self->flux = (double *)PyArray_DATA(self->df);
    self->ivar = (double *)PyArray_DATA(self->di);
    self->nobs = PyArray_DIMS(self->df)[0];
    self->nstars = PyArray_DIMS(self->df)[1];

    return self;
}

void PhotoData_destroy(PhotoData *self)
{
    if (self != NULL) {
        Py_XDECREF(self->df);
        Py_XDECREF(self->di);

        free(self);
    }
}

PhotoModel *PhotoModel_init(PyObject *model)
{
    PhotoModel *self = malloc(sizeof(PhotoModel));

    // parse the data/model classes
    self->mz = _fromNPYArray(model,"zero");
    self->mf = _fromNPYArray(model,"flux");
    if (self->mz == NULL || self->mf == NULL)
        return NULL;

    self->zero = (double *)PyArray_DATA(self->mz);
    self->flux = (double *)PyArray_DATA(self->mf);

    // model parameters
    self->jitterabs2 = _fromPyDouble(model,"jitterabs2");
    self->jitterrel2 = _fromPyDouble(model,"jitterrel2");
    self->Q2    = _fromPyDouble(model,"Q2");
    self->Pvar       = _fromPyDouble(model,"pvar");
    self->sigbad2    = _fromPyDouble(model,"sigbad2");
    self->Pbad       = _fromPyDouble(model,"pbad");

    return self;
}

void PhotoModel_destroy(PhotoModel *self)
{
    if (self != NULL) {
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
        PhotoData *data, PhotoModel *model)
{
    int ind = i*data->nstars+alpha; // data is TRANSPOSED from what we want!
    *lnpgood = 0.0; *lnpbad = 0.0; *lnpvargood = 0.0;

    if (data->ivar[ind] > 0.0) {
        double ff = model->zero[i]*model->flux[alpha]; // outer(zero,flux)
        double sig2   = 1.0/data->ivar[ind];
        double delta2 = model->jitterabs2 + model->jitterrel2*ff*ff;

        *lnpgood    = _lnnormal(data->flux[ind],ff,sig2+delta2);
        *lnpbad     = _lnnormal(data->flux[ind],ff,sig2+delta2+model->sigbad2);
        *lnpvargood = _lnnormal(data->flux[ind],ff,sig2+delta2+model->Q2*ff*ff);
    }
}

void get_lnpvar_and_lnpconst(int alpha, double *lnpvar, double *lnpconst,
        PhotoData *data, PhotoModel *model)
{
    int i;
    *lnpconst = 0.0;
    *lnpvar   = 0.0;
    for (i = 0; i < data->nobs; i++) {
        double lnpgood,lnpbad,lnpvargood;
        get_lnpgood_and_lnpbad_and_lnpvargood(i,alpha,
                &lnpgood,&lnpbad,&lnpvargood,
                data,model);
        *lnpconst += _logsumexp(log(1-model->Pbad)+lnpgood,
                                 log(model->Pbad)+lnpbad);
        *lnpvar   += _logsumexp(log(1-model->Pbad)+lnpvargood,
                                 log(model->Pbad)+lnpbad);
    }
}

double lnlikelihood(PhotoData *data, PhotoModel *model)
{
    int alpha;
    double lnlike = 0.0;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(alpha) \
    schedule(static) reduction(+:lnlike)
#endif
    for (alpha = 0; alpha < data->nstars; alpha++) {
        double lnpconst, lnpvar;
        get_lnpvar_and_lnpconst(alpha, &lnpvar, &lnpconst,
                data,model);
        lnlike += _logsumexp(log(1-model->Pvar)+lnpconst,log(model->Pvar)+lnpvar);
    }
    return lnlike;
}

void get_lnoddsvar(PhotoData *data, PhotoModel *model, double *lnoddsvar)
{
    int alpha;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(alpha) \
    schedule(static)
#endif
    for (alpha = 0; alpha < data->nstars; alpha++) {
        double lnpconst, lnpvar;

        get_lnpvar_and_lnpconst(alpha, &lnpvar, &lnpconst,
                data,model);

        lnoddsvar[alpha] = log(model->Pvar)-log(1-model->Pvar)
            + lnpvar - lnpconst;
    }
}

void get_lnlikeratiobad(PhotoData *data, PhotoModel *model, double *lnlikeratiobad)
{
    int i,alpha;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(alpha) \
    schedule(static)
#endif
    for (i = 0; i < data->nobs; i++) {
        for (alpha = 0; alpha < data->nstars; alpha++) {
            int ind = i*data->nstars+alpha;

            double lnpgood,lnpbad,lnpvargood;
            double lnpconst,lnpvar,lntotgood;

            get_lnpgood_and_lnpbad_and_lnpvargood(i,alpha,
                    &lnpgood,&lnpbad,&lnpvargood,
                    data,model);
            get_lnpvar_and_lnpconst(alpha, &lnpvar, &lnpconst,
                    data,model);

            //lntotgood = _logsumexp(lnpvar+lnpvargood,lnpconst+lnpgood);
            lntotgood = lnpgood;
            lnlikeratiobad[ind] = lnpbad - lntotgood;
        }
    }
}

void get_lnoddsbad(PhotoData *data, PhotoModel *model, double *lnoddsbad)
{
    int i,alpha;
#ifdef USEOPENMP
#pragma omp parallel for default(shared) private(alpha) \
    schedule(static)
#endif
    for (i = 0; i < data->nobs; i++) {
        for (alpha = 0; alpha < data->nstars; alpha++) {
            int ind = i*data->nstars+alpha;

            double lnpgood,lnpbad,lnpvargood;
            double lnpconst,lnpvar,lntotgood, lntotbad;

            get_lnpgood_and_lnpbad_and_lnpvargood(i,alpha,
                    &lnpgood,&lnpbad,&lnpvargood,
                    data,model);
            get_lnpvar_and_lnpconst(alpha, &lnpvar, &lnpconst,
                    data,model);

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
    PyObject *model0 = NULL, *data0 = NULL;
    PhotoData *data   = NULL;
    PhotoModel *model = NULL;
    if (!PyArg_ParseTuple(args, "OO", &model0, &data0))
        goto fail;

    data   = PhotoData_init(data0);
    model = PhotoModel_init(model0);
    if (data == NULL || model == NULL)
        goto fail;

    // Calcuate the likelihood
    double lnlike = lnlikelihood(data,model);
    PyObject *result = PyFloat_FromDouble(lnlike);

    // clean up!
    PhotoData_destroy(data);
    PhotoModel_destroy(model);
    return result;

fail:
    PyErr_SetString(PyExc_RuntimeError, "Incorrect input for likelihood");
    
    // clean up and fail
    PhotoData_destroy(data);
    PhotoModel_destroy(model);
    
    return NULL;
}

static PyObject *likelihood_lnoddsvar(PyObject *self, PyObject *args)
{
    PyObject *model0 = NULL, *data0 = NULL, *lnoddsvar_obj = NULL;
    PyObject *lnoddsvar = NULL;
    PhotoData *data   = NULL;
    PhotoModel *model = NULL;    
    if (!PyArg_ParseTuple(args, "OOO", &model0, &data0, &lnoddsvar_obj))
        goto fail;
    
    data   = PhotoData_init(data0);
    model = PhotoModel_init(model0);
    if (data == NULL || model == NULL)
        goto fail;

    lnoddsvar = PyArray_FROM_OTF(lnoddsvar_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    if (lnoddsvar == NULL)
        goto fail;
    get_lnoddsvar(data, model, (double *)PyArray_DATA(lnoddsvar));

    // clean up!
    Py_DECREF(lnoddsvar);
    PhotoData_destroy(data);
    PhotoModel_destroy(model);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    PyErr_SetString(PyExc_RuntimeError, "Incorrect input for lnoddsvar");
    
    Py_XDECREF(lnoddsvar);
    PhotoData_destroy(data);
    PhotoModel_destroy(model);
    
    return NULL;
}

static PyObject *likelihood_lnoddsbad(PyObject *self, PyObject *args)
{
    PyObject *model0 = NULL, *data0 = NULL, *lnoddsbad_obj = NULL, *lnoddsbad = NULL;
    PhotoData *data   = NULL;
    PhotoModel *model = NULL;
    if (!PyArg_ParseTuple(args, "OOO", &model0, &data0, &lnoddsbad_obj))
        goto fail;
    
    data = PhotoData_init(data0);
    model = PhotoModel_init(model0);
    if (data == NULL || model == NULL)
        goto fail;

    lnoddsbad = PyArray_FROM_OTF(lnoddsbad_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    if (lnoddsbad == NULL)
        goto fail;
    get_lnoddsbad(data, model, (double *)PyArray_DATA(lnoddsbad));

    // clean up!
    Py_DECREF(lnoddsbad);
    PhotoData_destroy(data);
    PhotoModel_destroy(model);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    PyErr_SetString(PyExc_RuntimeError, "Incorrect input for lnoddsbad");
    
    Py_XDECREF(lnoddsbad);
    PhotoData_destroy(data);
    PhotoModel_destroy(model);
    
    return NULL;
}

static PyObject *likelihood_lnlikeratiobad(PyObject *self, PyObject *args)
{
    PyObject *model0 = NULL, *data0 = NULL, *lnlikeratiobad_obj = NULL, *lnlikeratiobad = NULL;
    PhotoData *data   = NULL;
    PhotoModel *model = NULL;
    if (!PyArg_ParseTuple(args, "OOO", &model0, &data0, &lnlikeratiobad_obj))
        goto fail;
    
    data = PhotoData_init(data0);
    model = PhotoModel_init(model0);
    if (data == NULL || model == NULL)
        goto fail;

    lnlikeratiobad = PyArray_FROM_OTF(lnlikeratiobad_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    if (lnlikeratiobad == NULL)
        goto fail;
    get_lnlikeratiobad(data, model, (double *)PyArray_DATA(lnlikeratiobad));
    
    // clean up!
    Py_DECREF(lnlikeratiobad);
    PhotoData_destroy(data);
    PhotoModel_destroy(model);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    PyErr_SetString(PyExc_RuntimeError, "Incorrect input for lnlikeratiobad");
    
    Py_XDECREF(lnlikeratiobad);
    PhotoData_destroy(data);
    PhotoModel_destroy(model);
    
    return NULL;
}

