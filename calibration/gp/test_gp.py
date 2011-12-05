#!/usr/bin/env python
# encoding: utf-8
"""
Unit tests for the GP code

"""

import numpy as np
import gp

def test_grad():
    N = 50
    s = 10
    x = np.linspace(0,s,N)
    y = np.sin(x) + 0.1*np.random.randn(N)

    p = gp.GaussianProcess(a2=0.05, la2=0.05, s2=0.1)
    p.condition(x,y)
    grad0 = p.grad_loglike()
    grad1 = np.zeros(len(grad0))

    dx = 1e-4
    for i in range(len(grad0)):
        p[i] -= dx
        p.condition(x,y)
        lnlike1 = p.marginal_loglike()

        p[i] += 2*dx
        p.condition(x,y)
        lnlike2 = p.marginal_loglike()

        p[i] -= dx

        grad1[i] = (lnlike2 - lnlike1)/(2*dx)

    print grad0, grad1
    print grad0-grad1

if __name__ == '__main__':
    test_grad()

