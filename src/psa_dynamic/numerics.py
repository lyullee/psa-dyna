"""Sparse, bounded finite differences for stiff trace-species PSA states.

Unlike an adaptive numerical-Jacobian probe, these perturbations cannot grow
without bound for identically zero impurity states. Coloring is computed from
the declared dependency graph; no private SciPy functions are required.
"""
import numpy as np
from scipy.sparse import csc_matrix


def dense_gauss5_with_gauss3_error(function, accepted_times):
    """Integrate a dense trajectory per accepted solver interval.

    The five-point Gauss-Legendre estimate is paired with a three-point
    estimate. Summing absolute interval differences provides a conservative
    *screen*, not a rigorous bound for nonsmooth valve switches. The caller
    must fall back to adaptive quadrature if its requested tolerance fails.
    """
    times = np.asarray(accepted_times, dtype=float)
    if times.ndim != 1 or len(times) < 2 or not np.all(np.diff(times) > 0):
        raise ValueError("Accepted times must be strictly increasing")
    nodes5, weights5 = np.polynomial.legendre.leggauss(5)
    nodes3, weights3 = np.polynomial.legendre.leggauss(3)
    integral = None
    error_sum = None
    for start, end in zip(times[:-1], times[1:]):
        mid = 0.5 * (start + end)
        half = 0.5 * (end - start)
        values5 = [np.asarray(function(mid + half * node), dtype=float) for node in nodes5]
        values3 = [np.asarray(function(mid + half * node), dtype=float) for node in nodes3]
        if any(value.shape != values5[0].shape for value in values5 + values3):
            raise ValueError("Dense quadrature function changed output shape")
        five = half * np.tensordot(weights5, values5, axes=(0, 0))
        three = half * np.tensordot(weights3, values3, axes=(0, 0))
        if integral is None:
            integral = np.zeros_like(five)
            error_sum = np.zeros_like(five)
        integral += five
        error_sum += np.abs(five - three)
    if not np.isfinite(integral).all() or not np.isfinite(error_sum).all():
        raise FloatingPointError("Non-finite dense quadrature result")
    return integral, error_sum


class BoundedJacobian:
    def __init__(self, function, pattern, scales, relative_step=None):
        self.function=function
        self.pattern=csc_matrix(pattern)
        self.scales=np.asarray(scales,float)
        self.relative_step=np.sqrt(np.finfo(float).eps) if relative_step is None else float(relative_step)
        if not 1e-12<=self.relative_step<=1e-3:
            raise ValueError("Finite-difference relative step must lie between 1e-12 and 1e-3")
        if self.pattern.shape!=(len(self.scales),len(self.scales)) or np.any(self.scales<=0):
            raise ValueError("Jacobian pattern and strictly positive scales must agree")
        rows=[set(self.pattern.indices[self.pattern.indptr[i]:self.pattern.indptr[i+1]]) for i in range(len(self.scales))]
        self.groups=[];occupied=[]
        for col in sorted(range(len(rows)),key=lambda i:len(rows[i]),reverse=True):
            if not rows[col]:
                continue
            for group,used in zip(self.groups,occupied):
                if used.isdisjoint(rows[col]):
                    group.append(col);used.update(rows[col]);break
            else:
                self.groups.append([col]);occupied.append(set(rows[col]))

    def __call__(self,time,state):
        base=self.function(time,state)
        step=self.relative_step*np.maximum(abs(state),self.scales)
        data=np.zeros(self.pattern.nnz)
        for group in self.groups:
            trial=state.copy();trial[group]+=step[group]
            difference=self.function(time,trial)-base
            for col in group:
                start,stop=self.pattern.indptr[col:col+2]
                indices=self.pattern.indices[start:stop]
                data[start:stop]=difference[indices]/(trial[col]-state[col])
        return csc_matrix((data,self.pattern.indices,self.pattern.indptr),shape=self.pattern.shape)


class ComplexStepJacobian(BoundedJacobian):
    """Branch-local sparse derivatives with no subtraction of nearby RHSs.

    Real-valued branch selection in the physical core fixes the active valve/
    upwind branch. At a switching surface this is a one-branch derivative,
    not a claim that the original one-way law is globally differentiable.
    """
    def __call__(self,time,state):
        step=1e-20*self.scales
        data=np.zeros(self.pattern.nnz)
        for group in self.groups:
            trial=np.asarray(state,dtype=complex).copy()
            trial[group]+=1j*step[group]
            imaginary=np.imag(self.function(time,trial))
            for col in group:
                start,stop=self.pattern.indptr[col:col+2]
                indices=self.pattern.indices[start:stop]
                data[start:stop]=imaginary[indices]/step[col]
        return csc_matrix((data,self.pattern.indices,self.pattern.indptr),shape=self.pattern.shape)
