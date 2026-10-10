"""Explicit unit-bearing configuration for the frozen Phase 0 benchmark."""
from dataclasses import dataclass, asdict
import math


class NumericalError(RuntimeError):
    """An invalid or unresolved numerical result; never a feasible response."""
    def __init__(self, message, code='numerical_failure'):
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class Material:
    """Synthetic isotropic material: E [N/mm²], density [kg/mm³], t [mm]."""
    E: float = 70000.0
    nu: float = 0.30
    density: float = 2.70e-6
    thickness: float = 5.0

    def __post_init__(self):
        if not all(math.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Material values must be finite')
        if self.E <= 0 or self.density <= 0 or self.thickness <= 0:
            raise ValueError('E, density and thickness must be positive')
        if not -1 < self.nu < .5:
            raise ValueError('Stable isotropic material requires -1 < nu < 0.5')


@dataclass(frozen=True)
class Benchmark:
    """Defaults are PROMETHEUS-01-v0; overrides must be archived explicitly."""
    volume_fraction: float = .40
    load: float = 1000.0
    filter_radius: float = 6.0
    p: float = 3.0
    emin_ratio: float = 1e-6
    xmin: float = 1e-3
    move: float = .20
    max_evaluations: int = 300
    min_evaluations: int = 30

    def __post_init__(self):
        if not all(math.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Benchmark values must be finite')
        if not 0 < self.volume_fraction <= 1:
            raise ValueError('Volume fraction must lie in (0, 1]')
        if self.load <= 0 or self.filter_radius <= 0 or self.p < 1:
            raise ValueError('Positive load/radius and p >= 1 required')
        if not 0 < self.emin_ratio < 1 or not 0 < self.xmin < 1:
            raise ValueError('Void modulus ratio and xmin must lie in (0, 1)')
        if not 0 < self.move <= 1:
            raise ValueError('Move must lie in (0, 1]')
        if any(isinstance(v, bool) or not isinstance(v, int) or v < 1
               for v in (self.max_evaluations, self.min_evaluations)):
            raise ValueError('Evaluation limits must be positive integers')


@dataclass(frozen=True)
class LoadCase:
    """Magnitude [N]; angle [degrees] measured from DOWN toward RIGHT."""
    magnitude: float = 1000.0
    angle_degrees: float = 0.0

    def __post_init__(self):
        if not math.isfinite(self.magnitude) or self.magnitude <= 0:
            raise ValueError('Load magnitude must be positive and finite')
        if not math.isfinite(self.angle_degrees):
            raise ValueError('Load angle must be finite')

    def coefficients(self):
        theta = math.radians(self.angle_degrees)
        return (self.magnitude * math.sin(theta), self.magnitude * math.cos(theta))
