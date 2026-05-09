"""Binary wrapper"""

from __future__ import annotations

import logging
import subprocess
from abc import ABC, abstractmethod
from functools import wraps
from pathlib import Path
from typing import Callable, Literal, TypeVar

try:
    # Try the modern location first
    from typing import Concatenate, ParamSpec
except ImportError:
    # Fallback to the extension library
    from typing_extensions import Concatenate, ParamSpec


from ..lib import SeqTypes
from ..lib._utils import CheckAttrs

_T = TypeVar("_T")
_P = ParamSpec("_P")
_R = TypeVar("_R")


def _check_attributes(*attrs: str) -> Callable[[Callable[Concatenate[_T, _P], _R]], Callable[Concatenate[_T, _P], _R]]:
    """Decorator to ensure specific attributes are initialized before executing the function.

    Args:
        *attrs: Attribute names to check in the instance.

    Raises:
        AttributeError: If any specified attribute is `False` in the instance.
    """
    var_mapping = {"done": "run"}
    invalid_attrs = [attr for attr in attrs if attr not in var_mapping]
    if invalid_attrs:
        raise AttributeError(f"Invalid attribute names: {invalid_attrs}")

    def decorator(func: Callable[Concatenate[_T, _P], _R]) -> Callable[Concatenate[_T, _P], _R]:
        @wraps(func)
        def wrapper(instance: _T, *args: _P.args, **kwargs: _P.kwargs) -> _R:
            """Validate variable inequality and execute the wrapped function."""
            false_attrs = CheckAttrs.is_false(instance, *attrs)
            for var in sorted(false_attrs, key=lambda x: list(var_mapping.keys()).index(x)):
                raise AttributeError(f"Please run the {var_mapping[var]} method first.")
            return func(instance, *args, **kwargs)

        return wrapper

    return decorator


class BinaryWrapper(ABC):
    _prog: str
    _cmd_log: Literal["stdout", "stderr"] = "stdout"
    __slots__ = ("_logger", "_output", "_cmd", "_result", "done")

    def __init__(self, file: str | Path, output: str | Path | None = None, *args, **kwargs) -> None:
        file = Path(file)
        if not file.exists():
            raise FileNotFoundError(f"{file}")
        self._output = Path(output) if output else None
        args, kwargs = self._params_check(*args, **kwargs)
        self._construct_cmd(file, self._output, *args, **kwargs)
        self._cmd: list[str]
        self.done = False

    def __init_subclass__(cls, **kwargs) -> None:
        super().__init_subclass__(**kwargs)
        cls._logger = logging.getLogger(f"{cls.__module__}.{cls.__name__}")

    def run(self) -> None:
        """Execute the command."""
        if self._output:
            self._output.parent.mkdir(parents=True, exist_ok=True)
        self._logger.debug(self.cmd)
        try:
            result = subprocess.run(self._cmd, capture_output=True, check=True, text=True)
            self._logger.debug("%s", getattr(result, self._cmd_log))
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"{self._prog} failed with cmd: {self.cmd}\n{e.stderr}")

        self._result = self._output if self._output else result.stdout
        self._post_run()
        self.done = True

    @property
    @_check_attributes("done")
    def result(self) -> str | Path:
        return self._result

    @property
    def cmd(self) -> str:
        return " ".join(self._cmd)

    def _params_check(self, *args, **kwargs) -> tuple[tuple, dict]:
        return args, kwargs

    @abstractmethod
    def _construct_cmd(self, file: Path, output: Path | None, *args, **kwargs) -> None: ...

    def _post_run(self):
        pass


class TreeToolWrapper(BinaryWrapper):
    __slots__ = ("_model",)

    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *args,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        model: str = "AUTO",
        **kwargs,
    ) -> None:
        super().__init__(file, output, *args, seqtype=seqtype, model=model, **kwargs)
        self._model: str = model

    @property
    @_check_attributes("done")
    def model(self) -> str:
        return self._model

    def _params_check(self, *args, seqtype: str, **kwargs) -> tuple[tuple, dict]:
        if seqtype == SeqTypes.DNA:
            seqtype = "DNA"
        elif seqtype == SeqTypes.PEP:
            seqtype = "AA"
        else:
            seqtype = "AUTO"
        return super()._params_check(*args, seqtype=seqtype, **kwargs)

    @abstractmethod
    def _construct_cmd(
        self,
        file: Path,
        output: Path,
        *args,
        seqtype: Literal["DNA", "AA"] | None,
        model: str,
        **kwargs,
    ) -> None: ...
