"""Exact angular compliance envelope for a frozen linear two-load system.

Q = F.T @ solve(K,F), F columns RIGHT and DOWN, with the desired load
magnitude already included. C(theta) = a.T Q a, a=[sin(theta),cos(theta)].
Angles are radians. This is exact only for the supplied discretized linear
model and fixed load patch. It does not cover changing K/supports/contact.
The Phase 0 utility is retained unchanged as a historical calculation.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class EnvelopeResult:
    maximum: float
    angle_radians: float
    candidates: tuple


def compliance_envelope(q, low, high):
    q = np.asarray(q, dtype=float)
    if q.shape != (2, 2) or not np.all(np.isfinite(q)):
        raise ValueError('Q must be a finite 2 by 2 matrix')
    if not np.isfinite(low) or not np.isfinite(high) or low > high:
        raise ValueError('Angle interval must be finite and ordered')
    scale = max(float(np.linalg.norm(q)), np.finfo(float).tiny)
    if np.linalg.norm(q-q.T) > 1e-12 * scale:
        raise ValueError('Q must be symmetric; inconsistent solved loads')
    q = (q+q.T)/2  # only roundoff after the explicit symmetry check
    if np.linalg.eigvalsh(q)[0] < -1e-12 * scale:
        raise ValueError('Compliance Gram matrix must be positive semidefinite')
    # Period is pi, so an interval wider than pi contains a full period.
    search_high = min(float(high), float(low)+np.pi)
    a, b = q[0, 0]-q[1, 1], 2*q[0, 1]
    angles = [float(low), float(search_high)]
    if abs(a)+abs(b) > 0:
        base = np.arctan2(-b, a)
        first = int(np.floor((2*low-base)/np.pi))-1
        last = int(np.ceil((2*search_high-base)/np.pi))+1
        angles += [float((base+k*np.pi)/2) for k in range(first,last+1)
                   if low <= (base+k*np.pi)/2 <= search_high]
    values = [float(np.array([np.sin(t),np.cos(t)]) @ q @
                    np.array([np.sin(t),np.cos(t)])) for t in angles]
    best = int(np.argmax(values))
    return EnvelopeResult(values[best], angles[best], tuple(angles))
