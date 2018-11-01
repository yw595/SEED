#include <Python.h>
#include <modsupport.h>
#include <stdio.h>

/*
 * Function to be called from Python
 */
static PyObject* py_myFunction(PyObject* self, PyObject* args)
{
  char str[1000];
  char* x;
  char* y;
  PyArg_ParseTuple(args, "ss", &x, &y);
  strcpy(str,"/mnt/vdb/home/ubuntu2/miniconda2/bin/glpsol --lp ");
  strcat(str,x);
  strcat(str," --tmlim 5 -o ");
  strcat(str,y);
  int status = system(str);
  return Py_BuildValue("i", status);
}

/*
 * Another function to be called from Python
 */
static PyObject* py_myOtherFunction(PyObject* self, PyObject* args)
{
  char str[1000];
  char* x;
  char* y;
  char* z;
  PyArg_ParseTuple(args, "sss", &x, &y, &z);
  strcpy(str,"/mnt/vdb/home/ubuntu2/minDisj ");
  strcat(str,x);
  strcat(str," ");
  strcat(str,y);
  strcat(str," > ");
  strcat(str,z);
  int status = system(str);
  /* printf("%d", status); */
  /* fflush(stdout); */
  return Py_BuildValue("i",status);
}

/*
 * Bind Python function names to our C functions
 */
static PyMethodDef testpython_methods[] = {
  {"myFunction", py_myFunction, METH_VARARGS},
  {"myOtherFunction", py_myOtherFunction, METH_VARARGS},
  {NULL, NULL}
};

/*
 * Python calls this to let us initialize our module
 */
void inittestpython()
{
  (void) Py_InitModule("testpython", testpython_methods);
}
