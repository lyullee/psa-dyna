"""Scalar/batched agreement checks, including active complex derivatives."""
from copy import deepcopy
import numpy as np
import pytest
from psa_dynamic.cases import wang_case,s015_case
from psa_dynamic.core import Bed,equilibrium,Network
from psa_dynamic.schedules import build_schedule
from psa_dynamic.vectorized import equilibrium_batch


def test_batch_equilibrium_matches_scalar_and_complex_derivatives():
    rng=np.random.default_rng(717)
    for case in (wang_case(n_cells=4),s015_case(n_cells=4)):
        ns=len(case["species"])
        for material in case["bed"]["materials"]:
            p=rng.uniform(.001,3.,(7,ns))+1e-20j*rng.normal(size=(7,ns))
            t=rng.uniform(299.,320.,7)+1e-20j*rng.normal(size=7)
            batch=equilibrium_batch(p,t,material)
            for j in range(7):
                scalar=equilibrium(p[j],t[j],material)
                for group in range(3):
                    np.testing.assert_allclose(batch[group][j].real,scalar[group].real,rtol=5e-14,atol=1e-14)
                    np.testing.assert_allclose(batch[group][j].imag/1e-20,scalar[group].imag/1e-20,rtol=1e-12,atol=1e-12)


@pytest.mark.parametrize("energy_mode",["internal_energy","source_sensible"])
@pytest.mark.parametrize("cells",[4,13])
def test_vectorized_faces_and_balances_equal_reference_for_both_kinetics(energy_mode,cells):
    rng=np.random.default_rng(238)
    for case in (wang_case(n_cells=cells),s015_case(n_cells=cells)):
        case["bed"]["energy_mode"]=energy_mode
        bed=Bed(case["bed"],case["species"])
        state=bed.initial(3.,np.full(bed.ns,1/bed.ns),310.).reshape(bed.n,bed.nv)
        state[:,:2*bed.ns]*=rng.uniform(.8,1.2,(bed.n,2*bed.ns))
        state[:,-2:]+=rng.uniform(-3,3,(bed.n,2))
        for complex_input in (False,True):
            x=state.astype(complex) if complex_input else state.copy()
            if complex_input:
                x+=1e-20j*rng.normal(size=x.shape)*np.maximum(abs(x),1e-4)
            new=bed.internal_flux(x.ravel());old=bed.internal_flux_scalar(x.ravel())
            for actual,expected in zip(new,old):
                np.testing.assert_allclose(actual.real,expected.real,rtol=1e-13,atol=1e-10)
                if complex_input:
                    np.testing.assert_allclose(actual.imag/1e-20,expected.imag/1e-20,rtol=2e-12,atol=1e-9)
            # Boundary fluxes are deliberately nonzero; both closures must
            # include exactly the same incoming/outgoing enthalpy.
            new[0][0]=.1/bed.ns;new[0][-1]=.07/bed.ns
            new[1][0]=200.;new[1][-1]=140.
            actual=bed.rates(x.ravel(),*new)
            expected=bed.rates_scalar(x.ravel(),*new)
            for a,e in zip(actual,expected):
                np.testing.assert_allclose(np.real(a),np.real(e),rtol=3e-13,atol=2e-10)
                if complex_input:
                    np.testing.assert_allclose(np.imag(a)/1e-20,np.imag(e)/1e-20,rtol=3e-12,atol=3e-9)


def test_full_network_matches_scalar_reference_in_every_slot():
    for case in (wang_case(n_cells=4),s015_case("equalized_8step",4)):
        net=Network(case,build_schedule(case["schedule"]))
        reference=deepcopy(net)
        reference.bed.internal_flux=reference.bed.internal_flux_scalar
        reference.bed.rates=reference.bed.rates_scalar
        y=net.yfeed*.99+.01/net.ns
        state=np.r_[np.concatenate([net.bed.initial(1.2+b,y,307.+b) for b in range(net.nb)]),np.zeros(net.nledger)]
        state[:net.process_size].reshape(net.nb,net.bed.n,net.bed.nv)[:,::2,:net.ns]*=1.01
        for slot in net.schedule["slots"]:
            np.testing.assert_allclose(net.rhs(0,state,slot),reference.rhs(0,state,slot),rtol=3e-13,atol=2e-10)


def test_explicit_unit_thermal_factors_are_backward_compatible():
    """The new source-local thermal fields must not alter legacy Wang results."""
    base = wang_case(n_cells=4)
    explicit = deepcopy(base)
    for material in explicit["bed"]["materials"]:
        material["energy_capacity_factor"] = 1.0
        material["conduction_factor"] = 1.0
    original_net = Network(base, build_schedule(base["schedule"]))
    explicit_net = Network(explicit, build_schedule(explicit["schedule"]))
    rng = np.random.default_rng(918)
    state = original_net.initial()
    state = state + rng.normal(0.0, 1.0e-4, state.size)
    assert np.allclose(explicit_net.bed.energy_capacity_factor, 1.0)
    assert np.allclose(explicit_net.bed.conduction_factor, 1.0)
    for slot in original_net.schedule["slots"]:
        np.testing.assert_allclose(
            explicit_net.rhs(0.0, state, slot),
            original_net.rhs(0.0, state, slot),
            rtol=1e-13,
            atol=1e-10,
        )
