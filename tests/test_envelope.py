import numpy as np
import pytest
from prometheus.envelope import compliance_envelope
from prometheus.models import Material, Benchmark, LoadCase


@pytest.mark.parametrize('q,expected,angle', [
    (np.eye(2),1.0,None), (np.diag([2.,1.]),1+np.sin(np.deg2rad(20))**2,None),
    (np.diag([1.,2.]),2.0,0),
    (np.array([[2.,.3],[.3,3.]]),(5+np.sqrt(1.36))/2,None)])
def test_analytic_envelopes(q,expected,angle):
    result=compliance_envelope(q,*np.deg2rad([-20,20]))
    assert result.maximum == pytest.approx(expected, rel=1e-13)
    if angle is not None:
        assert result.angle_radians == pytest.approx(angle,abs=1e-13)


def test_direct_solve_and_dense_grid():
    rng=np.random.default_rng(1729)
    for _ in range(40):
        a=rng.normal(size=(12,12)); k=a.T@a+2*np.eye(12)
        f=rng.normal(size=(12,2)); q=f.T@np.linalg.solve(k,f)
        result=compliance_envelope(q,*np.deg2rad([-20,20]))
        angles=np.linspace(*np.deg2rad([-20,20]),10001)
        weights=np.array([np.sin(angles),np.cos(angles)])
        direct=f@weights[:,::100]
        actual=np.einsum('ij,ij->j',direct,np.linalg.solve(k,direct))
        reduced=np.einsum('ij,ij->j',weights,q@weights)
        np.testing.assert_allclose(actual,reduced[::100],rtol=1e-12,atol=1e-13)
        assert result.maximum >= reduced.max()-1e-12*result.maximum
        assert (result.maximum-reduced.max())/result.maximum < 1e-7


def test_full_period_and_single_angle():
    q=np.array([[2.,.3],[.3,3.]])
    assert compliance_envelope(q,-20,20).maximum == pytest.approx(np.linalg.eigvalsh(q)[-1])
    assert compliance_envelope(q,0,0).maximum == pytest.approx(3)


@pytest.mark.parametrize('q,lo,hi',[
    (np.eye(3),0,1),(np.array([[1,1],[0,1]]),0,1),
    (-np.eye(2),0,1),(np.full((2,2),np.nan),0,1),
    (np.eye(2),2,1),(np.eye(2),0,np.inf)])
def test_invalid_inputs(q,lo,hi):
    with pytest.raises(ValueError):compliance_envelope(q,lo,hi)


def test_units_and_invalid_configuration():
    assert Material().E == 70000
    assert Benchmark().max_evaluations == 300
    assert LoadCase().coefficients() == (0,1000)
    for kw in ({'E':-1},{'nu':.5},{'density':0},{'thickness':float('nan')}):
        with pytest.raises(ValueError):Material(**kw)
    with pytest.raises(ValueError):Benchmark(volume_fraction=0)
    with pytest.raises(ValueError):Benchmark(max_evaluations=3.5)
