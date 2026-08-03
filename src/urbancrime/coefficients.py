"""Construction and validation of scalar FEM coefficients."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Real
from typing import Any, Callable, Mapping

import numpy as np

from .config import CoefficientSpec, ModelConfig


AnalyticProfile = Callable[[np.ndarray], np.ndarray]
_PROFILES: dict[str, Callable[[Mapping[str, Any]], AnalyticProfile]] = {}


def register_profile(name: str) -> Callable[[Callable[[Mapping[str, Any]], AnalyticProfile]], Callable]:
    def decorator(factory: Callable[[Mapping[str, Any]], AnalyticProfile]) -> Callable:
        if name in _PROFILES:
            raise ValueError(f"coefficient profile {name!r} is already registered")
        _PROFILES[name] = factory
        return factory

    return decorator


@register_profile("sinusoidal_x")
def _sinusoidal_x(parameters: Mapping[str, Any]) -> AnalyticProfile:
    offset = float(parameters["offset"])
    amplitude = float(parameters["amplitude"])
    wavelength = float(parameters["wavelength"])
    phase = float(parameters["phase"])
    if wavelength <= 0:
        raise ValueError("sinusoidal_x wavelength must be positive")
    return lambda x: offset + amplitude * np.sin(np.pi * x[0] / wavelength + phase)


@register_profile("cosine_y")
def _cosine_y(parameters: Mapping[str, Any]) -> AnalyticProfile:
    offset = float(parameters["offset"])
    amplitude = float(parameters["amplitude"])
    wavelength = float(parameters["wavelength"])
    phase = float(parameters["phase"])
    if wavelength <= 0:
        raise ValueError("cosine_y wavelength must be positive")
    return lambda x: offset + amplitude * np.cos(np.pi * x[1] / wavelength + phase)


@register_profile("gaussian_mixture")
def _gaussian_mixture(parameters: Mapping[str, Any]) -> AnalyticProfile:
    background = float(parameters["background"])
    components = tuple(parameters["components"])
    parsed: list[tuple[float, float, float, float]] = []
    for component in components:
        amplitude = float(component["amplitude"])
        center_x, center_y = (float(v) for v in component["center"])
        width = float(component["width"])
        if width <= 0:
            raise ValueError("Gaussian width must be positive")
        parsed.append((amplitude, center_x, center_y, width))

    def evaluate(x: np.ndarray) -> np.ndarray:
        result = np.full(x.shape[1], background, dtype=float)
        for amplitude, center_x, center_y, width in parsed:
            radius2 = (x[0] - center_x) ** 2 + (x[1] - center_y) ** 2
            result += amplitude * np.exp(-radius2 / width)
        return result

    return evaluate


@register_profile("diagonal_highway")
def _diagonal_highway(parameters: Mapping[str, Any]) -> AnalyticProfile:
    background = float(parameters["background"])
    amplitude = float(parameters["amplitude"])
    sharpness = float(parameters["sharpness"])
    intercept = float(parameters["intercept"])
    x_min = parameters.get("x_min")
    x_max = parameters.get("x_max")
    y_min = parameters.get("y_min")
    y_max = parameters.get("y_max", parameters.get("y_cutoff"))
    if sharpness <= 0:
        raise ValueError("diagonal_highway sharpness must be positive")

    def evaluate(x: np.ndarray) -> np.ndarray:
        ridge = amplitude * np.exp(-sharpness * (x[0] + x[1] - intercept) ** 2)
        active = np.ones(x.shape[1], dtype=bool)
        if x_min is not None:
            active &= x[0] > float(x_min)
        if x_max is not None:
            active &= x[0] < float(x_max)
        if y_min is not None:
            active &= x[1] > float(y_min)
        if y_max is not None:
            active &= x[1] < float(y_max)
        ridge = np.where(active, ridge, 0.0)
        return background + ridge

    return evaluate


@register_profile("piecewise_linear_ridge")
def _piecewise_linear_ridge(parameters: Mapping[str, Any]) -> AnalyticProfile:
    """Gaussian ridge following two lines separated at ``switch_x``.

    Each line is represented by ``a*x + b*y + c = 0``.  The profile is used
    for the two-segment I-90 geometry in the M3AS Chicago experiment, while
    remaining reusable for other piecewise-linear corridors.
    """

    background = float(parameters["background"])
    amplitude = float(parameters["amplitude"])
    sharpness = float(parameters["sharpness"])
    switch_x = float(parameters["switch_x"])
    left = tuple(float(value) for value in parameters["left_line"])
    right = tuple(float(value) for value in parameters["right_line"])
    if sharpness <= 0:
        raise ValueError("piecewise_linear_ridge sharpness must be positive")
    if len(left) != 3 or len(right) != 3:
        raise ValueError("piecewise ridge lines must contain [a, b, c]")

    def evaluate(x: np.ndarray) -> np.ndarray:
        left_distance = left[0] * x[0] + left[1] * x[1] + left[2]
        right_distance = right[0] * x[0] + right[1] * x[1] + right[2]
        distance = np.where(x[0] < switch_x, left_distance, right_distance)
        return background + amplitude * np.exp(-sharpness * distance**2)

    return evaluate


@register_profile("layered_y")
def _layered_y(parameters: Mapping[str, Any]) -> AnalyticProfile:
    breakpoints = np.asarray(parameters["breakpoints"], dtype=float)
    values = np.asarray(parameters["values"], dtype=float)
    if breakpoints.ndim != 1 or values.ndim != 1:
        raise ValueError("layered_y breakpoints and values must be one-dimensional")
    if values.size != breakpoints.size + 1:
        raise ValueError("layered_y requires one more value than breakpoints")
    if np.any(np.diff(breakpoints) <= 0):
        raise ValueError("layered_y breakpoints must be strictly increasing")

    def evaluate(x: np.ndarray) -> np.ndarray:
        return values[np.searchsorted(breakpoints, x[1], side="right")]

    return evaluate


@register_profile("scaled_gaussian_mixture")
def _scaled_gaussian_mixture(parameters: Mapping[str, Any]) -> AnalyticProfile:
    background = float(parameters["background"])
    x_scale = float(parameters["x_scale"])
    y_scale = float(parameters["y_scale"])
    components = tuple(parameters["components"])
    if x_scale <= 0 or y_scale <= 0:
        raise ValueError("scaled Gaussian coordinate scales must be positive")

    def evaluate(x: np.ndarray) -> np.ndarray:
        scaled_x, scaled_y = x_scale * x[0], y_scale * x[1]
        result = np.full(x.shape[1], background, dtype=float)
        for component in components:
            center_x, center_y = (float(value) for value in component["center"])
            width = float(component["width"])
            if width <= 0:
                raise ValueError("scaled Gaussian width must be positive")
            radius2 = (scaled_x - center_x) ** 2 + (scaled_y - center_y) ** 2
            result += float(component["amplitude"]) * np.exp(-radius2 / width)
        return result

    return evaluate


@register_profile("trigonometric_product")
def _trigonometric_product(parameters: Mapping[str, Any]) -> AnalyticProfile:
    offset = float(parameters["offset"])
    amplitude = float(parameters["amplitude"])
    x_wavelength = float(parameters["x_wavelength"])
    y_wavelength = float(parameters["y_wavelength"])
    x_frequency = float(parameters["x_frequency"])
    y_frequency = float(parameters["y_frequency"])
    x_phase = float(parameters["x_phase"])
    y_phase = float(parameters["y_phase"])
    x_function = str(parameters["x_function"])
    y_function = str(parameters["y_function"])
    functions = {"sin": np.sin, "cos": np.cos}
    if x_wavelength <= 0 or y_wavelength <= 0:
        raise ValueError("trigonometric wavelengths must be positive")
    if x_function not in functions or y_function not in functions:
        raise ValueError("trigonometric functions must be sin or cos")

    def evaluate(x: np.ndarray) -> np.ndarray:
        x_wave = functions[x_function](
            x_frequency * np.pi * x[0] / x_wavelength + x_phase
        )
        y_wave = functions[y_function](
            y_frequency * np.pi * x[1] / y_wavelength + y_phase
        )
        return offset + amplitude * x_wave * y_wave

    return evaluate


def analytic_profile(spec: CoefficientSpec) -> AnalyticProfile:
    if spec.kind != "analytic" or spec.profile is None:
        raise ValueError("analytic_profile requires an analytic CoefficientSpec")
    try:
        factory = _PROFILES[spec.profile]
    except KeyError as error:
        choices = ", ".join(sorted(_PROFILES))
        raise ValueError(f"unknown analytic profile {spec.profile!r}; available: {choices}") from error
    return factory(spec.parameters)


def as_coefficient(value: Any, function_space: Any, mesh: Any) -> Any:
    """Normalize a scalar, callable, Constant, or Function for branch-free forms."""

    from dolfinx import fem
    from petsc4py import PETSc

    if isinstance(value, Real):
        return fem.Constant(mesh, PETSc.ScalarType(float(value)))
    if isinstance(value, (fem.Constant, fem.Function)):
        return value
    if callable(value):
        result = fem.Function(function_space)
        result.interpolate(value)
        result.x.scatter_forward()
        return result
    raise TypeError(f"unsupported coefficient input: {type(value).__name__}")


def build_coefficient(spec: CoefficientSpec, function_space: Any, mesh: Any) -> Any:
    if spec.kind == "constant":
        return as_coefficient(float(spec.value), function_space, mesh)
    if spec.kind == "analytic":
        return as_coefficient(analytic_profile(spec), function_space, mesh)
    return _load_file_coefficient(spec, function_space)


def _load_file_coefficient(spec: CoefficientSpec, function_space: Any) -> Any:
    """Load a nodal array or an interpolated point cloud from ``.npy/.npz``."""

    from dolfinx import fem
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

    if spec.path is None:
        raise ValueError("file coefficient requires a path")
    source = spec.path
    result = fem.Function(function_space, name=spec.field_name or source.stem)
    if source.suffix == ".npy":
        values = np.asarray(np.load(source), dtype=float).reshape(-1)
        if values.size != result.x.array.size:
            raise ValueError(
                f"{source} contains {values.size} values, expected {result.x.array.size} local dofs"
            )
        result.x.array[:] = values
    elif source.suffix == ".npz":
        with np.load(source) as archive:
            field_name = spec.field_name or "values"
            if "points" not in archive or field_name not in archive:
                raise ValueError(f"{source} must contain 'points' and {field_name!r}")
            points = np.asarray(archive["points"], dtype=float)
            values = np.asarray(archive[field_name], dtype=float).reshape(-1)
        if points.ndim != 2 or points.shape[0] != values.size or points.shape[1] < 2:
            raise ValueError("point-cloud coefficient requires points=(n,2) and values=(n,)")
        linear = LinearNDInterpolator(points[:, :2], values)
        nearest = NearestNDInterpolator(points[:, :2], values)

        def evaluate(x: np.ndarray) -> np.ndarray:
            query = x[:2].T
            output = np.asarray(linear(query), dtype=float)
            missing = ~np.isfinite(output)
            if np.any(missing):
                output[missing] = nearest(query[missing])
            return output

        result.interpolate(evaluate)
    else:
        raise ValueError("file coefficients currently support only .npy and .npz")
    result.x.scatter_forward()
    return result


@dataclass(frozen=True)
class BurglaryCoefficients:
    eta: Any
    ast: Any
    q: Any


def build_model_coefficients(config: ModelConfig, function_space: Any, mesh: Any) -> BurglaryCoefficients:
    result = BurglaryCoefficients(
        eta=build_coefficient(config.eta, function_space, mesh),
        ast=build_coefficient(config.ast, function_space, mesh),
        q=build_coefficient(config.q, function_space, mesh),
    )
    validate_coefficient(result.eta, "eta", lower=0.0, strict_lower=True)
    validate_coefficient(result.q, "q", lower=0.0)
    return result


def validate_coefficient(
    coefficient: Any,
    name: str,
    *,
    lower: float | None = None,
    upper: float | None = None,
    strict_lower: bool = False,
) -> None:
    """Validate Constant/Function values on the current MPI rank."""

    raw = getattr(coefficient, "value", None)
    if raw is None and hasattr(coefficient, "x"):
        raw = coefficient.x.array
    values = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"coefficient {name} contains non-finite values")
    if lower is not None:
        invalid = values <= lower if strict_lower else values < lower
        if np.any(invalid):
            operator = ">" if strict_lower else ">="
            raise ValueError(f"coefficient {name} must be {operator} {lower}")
    if upper is not None and np.any(values > upper):
        raise ValueError(f"coefficient {name} must be <= {upper}")
