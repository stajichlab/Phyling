"""Binary wrapper"""

from __future__ import annotations

import logging
import subprocess
from abc import ABC
from functools import wraps
from pathlib import Path
from typing import Callable, Generic, Literal, TypeVar, Union, cast

try:
    # Try the modern location first
    from typing import Concatenate, ParamSpec
except ImportError:
    # Fallback to the extension library
    from typing_extensions import Concatenate, ParamSpec


from ..exception import AlreadyExecutedError, SeqtypeError
from ..lib import SeqTypes
from ..lib._utils import CheckAttrs

_T = TypeVar("_T")
_P = ParamSpec("_P")
_R = TypeVar("_R")
_O = TypeVar("_O", Path, Union[Path, None])
_S = TypeVar("_S", bound=Literal["DNA", "AA", "AUTO"])


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


class BinaryWrapper(ABC, Generic[_O]):
    _prog: str
    _cmd_log: Literal["stdout", "stderr"] = "stdout"
    __slots__ = ("_logger", "_file", "_output", "_cmd", "_result", "done")
    _output: _O

    def __init__(self, file: str | Path, output: str | Path | None = None) -> None:
        self._file = Path(file)
        if not self._file.exists():
            raise FileNotFoundError(f"{self._file}")
        self._output = cast(_O, Path(output) if output else None)
        self.done = False
        self._cmd: list[str]

    def __init_subclass__(cls, **kwargs) -> None:
        super().__init_subclass__(**kwargs)
        cls._logger = logging.getLogger(f"{cls.__module__}.{cls.__name__}")

    def run(self) -> None:
        """Execute the command."""
        if self.done:
            raise AlreadyExecutedError(
                f"The {self._prog} instance for {self._file.name} has already been executed. "
                "To run again, you must instantiate a new tool wrapper."
            )
        if not hasattr(self, "_cmd"):
            RuntimeError(f"{self._prog} failed without cmd constructed")
        if isinstance(self._output, Path):
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

    def _post_run(self):
        pass


class TreeToolWrapper(BinaryWrapper[Path], Generic[_S]):
    __slots__ = ("_seqtype", "_model")

    _ALLOWED_SEQTYPES: tuple[str, ...] = ("DNA", "AA", "AUTO")
    _seqtype: _S

    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        model: str = "AUTO",
        **kwargs,
    ) -> None:
        super().__init__(file, output)
        if seqtype == SeqTypes.DNA:
            target = "DNA"
        elif seqtype == SeqTypes.PEP:
            target = "AA"
        else:
            target = "AUTO"

        if target not in self._ALLOWED_SEQTYPES:
            raise SeqtypeError(
                f"Invalid seqtype: {seqtype} for {self.__class__.__name__}. Allowed choices are: {list(self._ALLOWED_SEQTYPES)}"
            )

        self._seqtype = cast(_S, target)
        self._model: str = model

    @property
    @_check_attributes("done")
    def model(self) -> str:
        return self._model
