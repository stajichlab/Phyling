"""Contains internal utilities and decorators for use within the package.

These are not part of the public API and may change without notice.
"""

from __future__ import annotations

import logging
import re
import shutil
import time
from functools import wraps
from pathlib import Path
from typing import Any, Callable, TypeVar
from zlib import crc32

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import NCBICodonTable, NCBICodonTableDNA, NCBICodonTableRNA
from Bio.Seq import Seq

from .. import AVAIL_CPUS
from ..exception import BinaryNotFoundError, SeqtypeError
from . import SeqTypes

_C = TypeVar("Callable", bound=Callable[..., Any])
logger = logging.getLogger(__name__)


def get_file_checksum(file: str | Path) -> int:
    """Computes the CRC checksum of a file.

    Args:
        file (str | Path): Path to the file whose checksum needs to be calculated.

    Returns:
        int: The computed CRC checksum of the file.

    Raises:
        FileNotFoundError: If the specified file does not exist.
    """
    file = Path(file)
    if not file.exists():
        raise FileNotFoundError(f"{file} not found.")
    crc = 0
    with open(file, "rb") as f:
        crc += crc32(f.read())
    return crc


def is_gzip_file(file: str | Path) -> bool:
    """Checks whether a file is in gzip format.

    This function reads the first two bytes of the file and checks if they match the gzip magic number.

    Args:
        file (str | Path): The path to the file to check.

    Returns:
        bool: True if the file is in gzip format, False otherwise.
    """
    with open(file, "rb") as f:
        magic_number = f.read(2)
    return magic_number == b"\x1f\x8b"


def guess_seqtype(seq: str, ignore_failed=False) -> str | None:
    """Guess the sequence type (seqtype) for a given sequence string.

    This function attempts to identify the type of a given sequence string. It returns "dna" for DNA coding sequences, "rna" for
    RNA sequences, "pep" for peptide sequences, or None if the seqtype cannot be determined and `ignore_failed` is set to True.

    Args:
        seq (str): The sequence string to analyze.
        ignore_failed (bool): If True, suppress exceptions when the seqtype cannot be determined. Defaults to False.

    Raises:
        SeqtypeError: if failed to guess the seqtype.

    Returns:
        str | None: "dna" for DNA sequences, "rna" for RNA sequences, "pep" for peptide sequences, or None if the seqtype cannot
        be determined and `ignore_failed` is True.
    """
    chars = set()
    chars.update(seq)
    if "-" in chars:
        chars.remove("-")
    if chars.issubset(set(NCBICodonTableDNA.nucleotide_alphabet + "N")):
        return SeqTypes.DNA
    elif chars.issubset(set(NCBICodonTableRNA.nucleotide_alphabet + "N")):
        return SeqTypes.RNA
    elif chars.issubset(set(NCBICodonTable.protein_alphabet + "BXZ")):
        return SeqTypes.PEP
    else:
        if ignore_failed:
            return None
        raise SeqtypeError("Cannot determine seqtype. Aborted.")


def remove_files(*files: Path) -> None:
    """Removes files.

    Args:
        *files (Path): Files to remove.
    """
    if files:
        for file in files:
            file.unlink(missing_ok=True)


def remove_dirs(*dirs: Path) -> None:
    """Removes directories.

    Args:
        *dirs (Path): Directories to remove.
    """
    if dirs:
        for dir in dirs:
            shutil.rmtree(dir)


def check_threads(func: _C) -> _C:
    """Decorator to validate and adjust the 'threads' parameter passed to a function.

    Ensures the 'threads' parameter is a positive integer and does not exceed the number of available CPU cores. If 'threads'
    exceeds the available cores, it is adjusted to the maximum available.

    Args:
        func (Callable): The function to wrap.

    Returns:
        Callable: The wrapped function with validated 'threads' parameter.

    Raises:
        TypeError: If the 'threads' argument is not an integer.
        ValueError: If the 'threads' argument is not greater than 0.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        """Validate and adjust the 'threads' parameter before passing to a function."""
        if "threads" in kwargs:
            kwargs["threads"] = _check(kwargs["threads"])
        else:
            args = list(args)
            import inspect

            sig = inspect.signature(func)
            params = list(sig.parameters)
            if "threads" in params:
                threads_index = params.index("threads")
                if len(args) > threads_index:
                    args[threads_index] = _check(args[threads_index])

        return func(*args, **kwargs)

    def _check(threads: int) -> int:
        if not isinstance(threads, int):
            raise TypeError("Argument threads must be an integer.")
        if threads <= 0:
            raise ValueError("Argument threads must be an integer greater than 0.")
        if threads > 0 and threads > AVAIL_CPUS:
            logger.warning(
                "Argument threads=%s exceeds available CPU cores %s. Setting it to %s.",
                threads,
                AVAIL_CPUS,
                AVAIL_CPUS,
            )
            threads = AVAIL_CPUS
        return threads

    return wrapper


class Timer:
    """Utility class for measuring elapsed time.

    Provides methods for starting and stopping a timer, as well as a decorator for measuring the execution time of a function.
    """

    __slots__ = ("_t_start",)

    def __init__(self) -> None:
        """Initialize the Timer instance."""
        self._t_start = None

    def start(self) -> None:
        """Start the timer.

        Raises:
            RuntimeError: If the timer is already running.
        """
        if self._t_start is not None:
            raise RuntimeError("Timer is already running.")
        self._t_start = time.monotonic()

    def stop(self) -> None:
        """Stop the timer and log the elapsed time as a formatted string.

        Raises:
            RuntimeError: If the timer has not been started.
        """
        if self._t_start is None:
            raise RuntimeError("Timer has not been started.")
        elapsed_time = time.monotonic() - self._t_start
        self._t_start = None
        logger.debug("Elapsed time: %s", self._format_elapsed(elapsed_time))

    @staticmethod
    def _format_elapsed(elapsed_sec: float) -> str:
        """Format the elapsed time into a human-readable string.

        Args:
            elapsed_sec (float): Elapsed time in seconds.

        Returns:
            str: Formatted string representing the elapsed time in hours, minutes, and seconds.
        """
        hours, remainder = divmod(elapsed_sec, 3600)
        mins, secs = divmod(remainder, 60)
        msg = f"{secs:.3f} seconds"
        if mins > 0:
            msg = f"{int(mins)} minutes " + msg
        if hours > 0:
            msg = f"{int(hours)} hours " + msg
        return msg

    @staticmethod
    def timer(func: _C) -> _C:
        """Decorator to measure the elapsed time of a function.

        Args:
            func (Callable): The function to be measured.

        Returns:
            Callable: The wrapped function that logs its execution time.
        """

        @wraps(func)
        def wrapper(*args, **kwargs):
            """Measure the elapsed time of a function and output to the logger."""
            t_start = time.monotonic()
            result = func(*args, **kwargs)
            elapsed_time = time.monotonic() - t_start
            logger.debug(
                "Elapsed time for %s: %s",
                func.__name__,
                Timer._format_elapsed(elapsed_time),
            )
            return result

        return wrapper


class CheckAttrs:
    """A utility class providing decorator-based checks for instance attributes.

    This class cannot be instantiated directly. Instead, use one of its nested classes (Exists, IsNotNone) as decorators to
    validate the presence or values of attributes in an object instance.
    """

    def __init__(self) -> None:
        raise TypeError("This class is not meant to be instantiated directly. Use it as a decorator.")

    class Exists:
        """A decorator class to ensure specified attributes exist on an instance.

        This class is used as a decorator to verify that certain attributes are present on an instance before executing a
        function. If any required attribute is missing, it raises an `AttributeError`.

        Attributes:
            attrs (tuple): A tuple of attribute names to check for existence.
        """

        __slots__ = ("attrs",)

        def __init__(self, *attrs: str) -> None:
            """A decorator to ensure the attributes exist on an instance.

            Args:
                *attrs (str): Attribute names to check for existence.
            """
            self.attrs = attrs

        def __call__(self, func: _C) -> _C:
            """Wrap the target function to perform attribute existence checks before execution.

            Args:
                func (Callable): The function to wrap.

            Returns:
                Callable: The wrapped function.

            Raises:
                AttributeError: If any required attribute is missing.
            """

            @wraps(func)
            def wrapper(instance, *args, **kwargs):
                """Validate attribute existence and execute the wrapped function."""
                missing_attrs = CheckAttrs.not_exists(instance, *self.attrs)
                if missing_attrs:
                    raise AttributeError(
                        "The following attributes are missing in the instance of "
                        f"{type(instance).__name__}: {', '.join(missing_attrs)}"
                    )
                return func(instance, *args, **kwargs)

            return wrapper

    class IsNotNone:
        """A decorator to ensure specified attributes of an instance are not None.

        This class wraps a function and ensures that the specified attributes of an instance are not `None` before executing the
        function. If any of the attributes are `None`, an `AttributeError` is raised.

        Attributes:
            attrs (tuple): A tuple of attribute names to check for None.
        """

        __slots__ = ("attrs",)

        def __init__(self, *attrs: str) -> None:
            """A decorator to ensure the specified attributes are not None on an instance.

            Args:
                *attrs (str): Attribute names to check for None.
            """
            self.attrs = attrs

        def __call__(self, func: _C) -> _C:
            """Wrap the target function to perform inequality checks before execution.

            Args:
                func (Callable): The function to wrap.

            Returns:
                Callable: The wrapped function.

            Raises:
                AttributeError: If any required attribute is None.
            """

            @wraps(func)
            def wrapper(instance, *args, **kwargs):
                """Validate attribute inequality and execute the wrapped function."""
                none_attrs = CheckAttrs.is_none(instance, *self.attrs)
                if none_attrs:
                    raise AttributeError(
                        f"The following attributes are None in the instance of {type(instance).__name__}: {', '.join(none_attrs)}"
                    )
                return func(instance, *args, **kwargs)

            return wrapper

    @staticmethod
    def not_exists(instance, *attrs: str) -> list[str]:
        """Check if all attributes exist on the given instance.

        Args:
            instance: The object to validate.
            *attrs (str): Attribute names to check.

        Returns:
            The attribute names that are missing from the instance.
        """
        missing_attrs = [attr for attr in attrs if not hasattr(instance, attr)]
        return missing_attrs

    @staticmethod
    def is_none(instance, *attrs: str) -> list[str]:
        """Check if any attributes are None on the given instance.

        Args:
            instance: The object to validate.
            *attrs (str): Attribute names to check.

        Returns:
            The attribute names that are None in the instance.
        """
        none_attrs = [attr for attr in attrs if getattr(instance, attr) is None]
        return none_attrs

    @staticmethod
    def is_false(instance, *attrs: str) -> list[str]:
        """Check if any attributes are False on the given instance.

        Args:
            instance: The object to validate.
            *attrs (str): Attribute names to check.

        Returns:
            The attribute names that are False in the instance.
        """
        none_attrs = [attr for attr in attrs if getattr(instance, attr) is False]
        return none_attrs


def check_binary(prog: str, bins: tuple, conda_url: str | None = None, source_url: str | None = None) -> str:
    """Find the path of a binary from the given sequence of entry points.

    This function searches for the binaries in the system's PATH and returns the path of the first binary found. If none are
    found, it raises an error.

    Args:
        *bins (str): Entry points (binary names) to search for.

    Returns:
        Path: The full path to the first available binary.

    Raises:
        BinaryNotFoundError: If none of the binaries are found in the system's PATH.
    """
    for bin in bins:
        bin_path = shutil.which(bin)
        if bin_path:
            return bin_path

    conda_msg = ""
    source_msg = ""
    if conda_url or source_url:
        install_msg = " Please install it "
    if conda_url:
        conda_msg = f"through 'conda install {conda_url}'"
        if source_url:
            conda_msg += " or "
    if source_url:
        source_msg = f"from source following the instruction on {source_url}"

    raise BinaryNotFoundError(f"{prog} not found.{install_msg}{conda_msg}{source_msg}")


def substitute_ambiguous_seq(seq: Seq) -> Seq:
    """Substitutes ambiguous characters in a sequence with a hyphen ('-') and converts the sequence to uppercase.

    Args:
        seq (Seq): The input sequence.

    Returns:
        Seq: The modified sequence with ambiguous characters replaced and converted to uppercase.
    """
    return Seq(re.sub(r"[ZzBbXx\*\.]", "-", str(seq))).upper()


def load_msa(file: Path) -> MultipleSeqAlignment:
    """Load a multiple sequence alignment (MSA) from a file in FASTA format.

    This function reads the MSA saved in FASTA format, processes each sequence to substitute ambiguous characters, and returns the
    alignment.

    Args:
        file (str | Path): Path to the FASTA file containing the multiple sequence alignment.

    Returns:
        MultipleSeqAlignment: The processed multiple sequence alignment object.

    Raises:
        ValueError: If the FASTA file is not in a valid format or cannot be read.
    """
    try:
        alignment = AlignIO.read(file, "fasta")
    except ValueError as e:
        raise ValueError(f"{file.absolute()}: {e}")
    for seq in alignment:
        seq.seq = substitute_ambiguous_seq(seq.seq)
    return alignment
