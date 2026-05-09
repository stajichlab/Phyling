"""Abstract classes for phyling."""

from __future__ import annotations

import textwrap
from abc import ABC, abstractmethod
from functools import wraps
from itertools import zip_longest
from pathlib import Path
from typing import Callable, Generic, Iterator, Literal, MutableSequence, Protocol, Sequence, TypeVar, overload

try:
    # Try the modern location first
    from typing import Concatenate, ParamSpec, Self
except ImportError:
    # Fallback to the extension library
    from typing_extensions import Concatenate, ParamSpec, Self

from ..exception import SeqtypeError
from . import SeqTypes
from ._utils import CheckAttrs, get_file_checksum

__all__ = [
    "FileWrapperABC",
    "SeqFileWrapperABC",
    "DataListABC",
    "SeqDataListABC",
]


_T = TypeVar("_T")
_P = ParamSpec("_P")
_R = TypeVar("_R")
_FileWrapperABC = TypeVar("_FileWrapperABC", bound="FileWrapperABC")
_SeqFileWrapperABC = TypeVar("_SeqFileWrapperABC", bound="SeqFileWrapperABC")
_DataListABC = TypeVar("_DataListABC", bound="DataListABC")


class SupportsSeqType(Protocol):
    """Structural type for any object with a seqtype attribute."""

    @property
    def seqtype(self) -> Literal["dna", "rna", "pep", "NaN"]: ...


_SupportsSeqType = TypeVar("_SupportsSeqType", bound=SupportsSeqType)


def check_class(func: Callable[_P, _R]) -> Callable[_P, _R]:
    """A decorator to ensure that the `other` argument is an instance of the same class as `instance`.

    Args:
        func (Callable): The function to be wrapped, which must take `instance` and `other` as arguments.

    Returns:
        Callable: A wrapped function that raises a `TypeError` if `other` is not of the same type as `instance`.

    Raises:
        TypeError: If `other` is not an instance of the same class as `instance`.

    Examples:
        @check_class
        def example_method(self, other):
            ...
    """

    @wraps(func)
    def wrapper(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        """Validate both objects are the same type."""
        if len(args) >= 2:
            if not isinstance(args[1], type(args[0])):
                raise TypeError(
                    f"Cannot operates with a different type. Expect a {type(args[1])} object but got a {type(args[0])}."
                )
        return func(*args, **kwargs)

    return wrapper


def check_loaded(
    func: Callable[Concatenate[_FileWrapperABC, _P], _R],
) -> Callable[Concatenate[_FileWrapperABC, _P], _R]:
    """A decorator to ensure that the data has been loaded before calling the decorated method.

    Args:
        func (Callable): The method to be wrapped, which operates on an instance with a `_data` attribute.

    Returns:
        Callable: A wrapped function that raises a `RuntimeError` if the `_data` attribute of the instance is `None`.

    Raises:
        RuntimeError: If the `_data` attribute of the instance is `None`, indicating that the required data has not been loaded.

    Examples:
        @check_loaded
        def process_data(self):
            ...
    """

    @wraps(func)
    def wrapper(instance: _FileWrapperABC, *args: _P.args, **kwargs: _P.kwargs) -> _R:
        """Validate whether the data has been loaded before executing the function."""
        if CheckAttrs.is_none(instance, "_data"):
            raise RuntimeError("Data is not loaded yet. Please run the load method to load it first.")
        return func(instance, *args, **kwargs)

    return wrapper


def load_data(
    func: Callable[Concatenate[_FileWrapperABC | _DataListABC, _P], _R],
) -> Callable[Concatenate[_FileWrapperABC | _DataListABC, _P], _R]:
    """A decorator to ensure that data is loaded before executing a method and unloaded afterward if it was not loaded initially.

    This decorator checks whether the `instance` object has its data loaded by evaluating `instance._data`. If the data is already
    loaded, the function executes without calling `load()` or `unload()`. If the data is not loaded, it calls the `load()` method
    before executing the wrapped function and the `unload()` method afterward.

    Args:
        func (Callable): The function to be wrapped, which requires the data to be loaded.

    Returns:
        Callable: A wrapped function that manages data loading and unloading.

    Raises:
        AttributeError: If `instance` does not have `load()` or `unload()` methods.

    Examples:
        @load_data
        def process_data(self, *args, **kwargs):
            ...
    """

    @wraps(func)
    def wrapper(instance: _FileWrapperABC | _DataListABC, *args: _P.args, **kwargs: _P.kwargs) -> _R:
        """Validate the data is loaded before executing the function and unload it afterward if it was not loaded initially."""
        was_loaded = instance._data
        if not was_loaded:
            instance.load()
        try:
            r = func(instance, *args, **kwargs)
        finally:
            if not was_loaded:
                instance.unload()
        return r

    return wrapper


def check_seqtype(instance: _SupportsSeqType, other: _SupportsSeqType) -> None:
    """Validate both objects represent the same sequence type."""
    if not other.seqtype == instance.seqtype:
        raise SeqtypeError("Items represent different seqtypes.")


def _add_seqtypes(x: _SupportsSeqType, y: _SupportsSeqType) -> Literal["dna", "rna", "pep", "NaN"]:
    """Determine the seqtype for the combined items x and y.

    Args:
        x (SeqFileWrapper | SeqDataList): The first item.
        y (SeqFileWrapper | SeqDataList): The second item.

    Returns:
        The seqtype if valid, or "NaN" if either item's seqtype is "NaN".

    Raises:
        SeqtypeError: If x and y have conflicting seqtypes.
    """
    seqtypes = {x.seqtype, y.seqtype}

    if len(seqtypes) == 1:
        return seqtypes.pop()

    if "NaN" in seqtypes:
        return (seqtypes - {"NaN"}).pop()

    raise SeqtypeError("Items represent different seqtypes.")


def list_repr_wrapper(data: list[_T]) -> str:
    """Format an iterator of data into a concise, multiline string representation.

    This function converts each element of the iterator to a string and limits the output to the first 5 and last 5 elements if
    the iterator has more than 20 elements. Excess elements are represented with an ellipsis ("..."). Each line is indented for
    readability.

    Args:
        data (Iterator): The data to be formatted.

    Returns:
        str: A string representation of the data, formatted for readability.

    Example:
        >>> data = range(25)
        >>> print(list_repr_wrapper(data))
        [0,
         1,
         2,
         3,
         4,
         ...,
         20,
         21,
         22,
         23,
         24]
    """
    str_data = [str(d) for d in data]
    if len(str_data) > 20:
        str_data = str_data[:5] + ["..."] + str_data[-5:]
    str_data = ",\n".join([textwrap.indent(line, prefix=" ") for line in str_data]).lstrip(" ")
    return f"[{str_data}]"


class FileWrapperABC(ABC):
    """An abstract base class for representing and managing files.

    This class provides a framework for handling file objects with attributes like file path, name, and checksum, as well as
    methods for file comparison and management.

    Attributes:
        file (pathlib.Path): The path to the file.
        name (str): The representative name of the file.
        checksum (int): The CRC checksum of the file, used for integrity verification.
    """

    __slots__ = ("_file", "_checksum", "_name", "_data")

    def __init__(self, file: str | Path, name: str | None = None) -> None:
        """Initialize the object with a file path and an representative name.

        Args:
            file (str | Path): The path to the file. It will be converted to a `pathlib.Path` object.
            name (str | None, optional): The representative name of the file. Defaults to the file name.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            RuntimeError: If the specified path is not a file.
        """
        self.file = file
        self.name = name
        self._data: MutableSequence = []

    def __repr__(self) -> str:
        """Return a string representation of the object.

        Returns:
            str: The string representation in the format `<ClassName(file_name)>`.
        """
        return f"{type(self).__qualname__}({self.name})"

    @check_class
    def __gt__(self, other: Self) -> bool:
        """Compare if the current object is larger (posterior) than another object.

        Args:
            other (FileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the current object's name is larger, otherwise False.
        """
        return self.name > other.name

    @check_class
    def __ge__(self, other: Self) -> bool:
        """Compare if the current object is larger than or equal to another object.

        Args:
            other (FileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the current object's name is larger or equal, otherwise False.
        """
        return self.name >= other.name

    @check_class
    def __lt__(self, other: Self) -> bool:
        """Compare if the current object is smaller (prior) than another object.

        Args:
            other (FileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the current object's name is smaller, otherwise False.
        """
        return self.name < other.name

    @check_class
    def __le__(self, other: Self) -> bool:
        """Compare if the current object is smaller than or equal to another object.

        Args:
            other (FileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the current object's name is smaller or equal, otherwise False.
        """
        return self.name <= other.name

    def __eq__(self, other: object) -> bool:
        """Check if the current object is equal to another object.

        Args:
            other (FileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the objects have the same name, otherwise False.
        """
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return self.checksum

    @property
    def file(self) -> Path:
        """Get the file path.

        Returns:
            Path: The file path as a `pathlib.Path` object.
        """
        return self._file

    @file.setter
    def file(self, file) -> None:
        """Set the file path."""
        file = Path(file).absolute()
        if not file.exists():
            raise FileNotFoundError(f"{self._file}")
        if not file.is_file():
            raise RuntimeError(f"{self._file} is not a file.")
        self._file = file
        self._checksum = get_file_checksum(self._file)

    @property
    def name(self) -> str:
        """Get the representative name of the file.

        Returns:
            str: The name of the file.
        """
        return self._name

    @name.setter
    def name(self, name) -> None:
        """Set the representative name of the file."""
        self._name = name if name else self.file.name

    @property
    def checksum(self) -> int:
        """Get the CRC checksum of the file.

        Returns:
            int: The checksum value of the file.
        """
        return self._checksum

    @abstractmethod
    def load(self) -> None:
        """Load the file content.

        This method must be implemented by subclasses to handle the specific logic for loading files.
        """
        ...

    def unload(self) -> None:
        """Clear the loaded content from memory.

        This resets the `_data` attribute to empty list.
        """
        self._data = []


class SeqFileWrapperABC(FileWrapperABC):
    """An abstract base class for representing and managing sequence files.

    This class provides a framework for handling sequence file objects, including attributes such as file path, representative
    name, CRC checksum, and sequence type. It also offers methods for file comparison, validation, and management.

    Attributes:
        file (pathlib.Path): The path to the sequence file.
        name (str): The representative name of the file.
        checksum (int): The CRC checksum of the file, used for integrity verification.
        seqtype (str): The sequence type of the file (e.g., DNA, RNA, Protein), determined automatically.
    """

    __slots__ = ("_seqtype",)
    __hash__ = FileWrapperABC.__hash__

    def __init__(self, file: str | Path, name: str | None = None, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None:
        """Initialize the object with a file path and an representative name.

        Args:
            file (str | Path): The path to the file. It will be converted to a `pathlib.Path` object.
            name (str | None, optional): The representative name of the file. Defaults to the file name.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            RuntimeError: If the specified path is not a file.
        """
        super().__init__(file, name)
        if seqtype == "AUTO":
            self._seqtype = self._guess_seqtype()
        elif seqtype == SeqTypes.DNA:
            self._seqtype = SeqTypes.DNA
        elif seqtype == SeqTypes.PEP:
            self._seqtype = SeqTypes.PEP
        else:
            raise SeqtypeError(f"Invalid seqtype: {seqtype}.")

    def __repr__(self) -> str:
        """Return a string representation of the object.

        Returns:
            str: The string representation in the format `<ClassName(file_name; seqtype=seqtype)>`.
        """
        return super().__repr__()[:-1] + f"; seqtype={self.seqtype})"

    def __gt__(self, other: Self) -> bool:
        """Compare if the current object is larger (posterior) than another object.

        Args:
            other (SeqFileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the current object's name is larger, otherwise False.

        Raises:
            SeqtypeError: If `instance.seqtype` and `other.seqtype` are not identical.
        """
        check_seqtype(self, other)
        return super().__gt__(other)

    def __ge__(self, other: Self) -> bool:
        """Compare if the current object is larger than or equal to another object.

        Args:
            other (SeqFileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the current object's name is larger or equal, otherwise False.

        Raises:
            SeqtypeError: If `instance.seqtype` and `other.seqtype` are not identical.
        """
        check_seqtype(self, other)
        return super().__ge__(other)

    def __lt__(self, other: Self) -> bool:
        """Compare if the current object is smaller (prior) than another object.

        Args:
            other (SeqFileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the current object's name is smaller, otherwise False.

        Raises:
            SeqtypeError: If `instance.seqtype` and `other.seqtype` are not identical.
        """
        check_seqtype(self, other)
        return super().__lt__(other)

    def __le__(self, other: Self) -> bool:
        """Compare if the current object is smaller than or equal to another object.

        Args:
            other (SeqFileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the current object's name is smaller or equal, otherwise False.

        Raises:
            SeqtypeError: If `instance.seqtype` and `other.seqtype` are not identical.
        """
        check_seqtype(self, other)
        return super().__le__(other)

    def __eq__(self, other: object) -> bool:
        """Check if the current object is equal to another object.

        Args:
            other (SeqFileWrapperABC): Another instance of the same class.

        Returns:
            bool: True if the objects have the same name, otherwise False.

        Raises:
            SeqtypeError: If `instance.seqtype` and `other.seqtype` are not identical.
        """
        if not isinstance(other, type(self)):
            return NotImplemented
        check_seqtype(self, other)
        return super().__eq__(other)

    @property
    def seqtype(self) -> Literal["dna", "pep", "rna", "NaN"]:
        """Get the sequence type of the file.

        Returns:
            str: The sequence type of the file.
        """
        return self._seqtype

    @abstractmethod
    def _guess_seqtype(self) -> Literal["dna", "pep", "rna", "NaN"]:
        """Guess the sequence type of the file.

        This method must be implemented by subclasses to handle the specific logic for guessing sequence type.
        """
        ...


class DataListABC(ABC, Generic[_FileWrapperABC]):
    """A list-like abstract base class providing partial list methods for managing FileWrapperABC objects and their associated
    metadata.

    Attributes:
        files (tuple): A tuple of file paths.
        names (tuple): A tuple of file names.
        checksums (dict): A dictionary of names and their checksums.
    """

    __slots__ = ("_data",)
    _bound_class: type[_FileWrapperABC]

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, data: Sequence[str | Path | _FileWrapperABC]) -> None: ...
    @overload
    def __init__(self, data: Sequence[str | Path | _FileWrapperABC], names: Sequence[str], *args, **kwargs) -> None: ...
    def __init__(self, data: Sequence[str | Path | _FileWrapperABC] = (), names: Sequence[str] = (), *args, **kwargs) -> None:
        """Initializes the object and stores data into a list.

        Args:
            data (Sequence[str | Path | FileWrapperABC], optional): A sequence of data items.
            names (Sequence[str], optional): A sequence of names corresponding to the data items.

        Raises:
            RuntimeError: If names are provided but data is not.
            TypeError: If a data item cannot be converted to the bound class.
            KeyError: If the item already exists.
        """
        self._data: list[_FileWrapperABC] = []

        if names and not data:
            raise RuntimeError("Received no data with names specified.")

        if data:
            if names and len(data) != len(names):
                raise RuntimeError("Data and names have different length.")

            for d, name in zip_longest(data, names or []):
                if isinstance(d, (str, Path)):
                    d = self._bound_class(d, name, *args, **kwargs)
                if not isinstance(d, self._bound_class):
                    raise TypeError(f"{type(d).__qualname__} cannot be converted to {self._bound_class.__qualname__}.")
                self.append(d)
        else:
            if names:
                raise RuntimeError("Received no data with names specified.")

    def __repr__(self) -> str:
        """Returns a string representation of the object.

        Returns:
            str: The string representation of the object.
        """
        return type(self).__qualname__ + "\n" + list_repr_wrapper(self._data)

    def __len__(self) -> int:
        """Returns the number of items in the object.

        Returns:
            int: The number of items in the object.
        """
        return len(self._data)

    def __iter__(self) -> Iterator[_FileWrapperABC]:
        """Returns an iterator over the items in the object.

        Returns:
            Iterator[FileWrapperABC]: An iterator for the items.
        """
        return iter(self._data)

    @check_class
    def __eq__(self, other: object) -> bool:
        """Checks if two objects have the same set of names.

        Args:
            other (DataListABC): Another object to compare.

        Returns:
            bool: True if the sets of names are equal, False otherwise.

        Raises:
            TypeError: If `other` is not an instance of the same class as `self`.
        """
        if not isinstance(other, type(self)):
            return NotImplemented

        return set(self.names) == set(other.names)

    def __contains__(self, item: str | _FileWrapperABC) -> bool:
        """Checks if an item is present in the object.

        Args:
            item (FileWrapperABC): The item to check.

        Returns:
            bool: True if the item is present, False otherwise.

        Raises:
            TypeError: If the item is not a string or a file wrapper object.
        """
        if isinstance(item, str):
            name = item
        elif isinstance(item, FileWrapperABC):
            name = item.name
        else:
            raise TypeError(f"Can only check by str or {FileWrapperABC.__qualname__} and its subclass.")
        return name in self.names

    @overload
    def __getitem__(self, key: str) -> _FileWrapperABC: ...
    @overload
    def __getitem__(self, key: int) -> _FileWrapperABC: ...
    @overload
    def __getitem__(self, key: slice) -> Self: ...
    def __getitem__(self, key: str | int | slice) -> _FileWrapperABC | Self:
        """Retrieves an item or subset of items by name, index, or slice.

        Args:
            key (str | int | slice): The key to retrieve.

        Returns:
            FileWrapperABC | DataListABC: The corresponding item or subset of items.
        """
        if isinstance(key, str):
            if key not in self.names:
                raise KeyError(f"{key}: Sample not found.")
            return self._data[self.names.index(key)]
        elif isinstance(key, slice):
            return self.__class__(self._data[key])
        else:
            return self._data[key]

    def __add__(self, other: Self) -> Self:
        """Concatenates two DataListABC objects.

        Args:
            other (DataListABC): The other DataListABC object.

        Returns:
            DataListABC: A new DataListABC object containing the concatenated data.
        """
        return self.__class__(self._data + other._data, self.names + other.names)

    @property
    def files(self) -> tuple[Path, ...]:
        """Returns the file paths in pathlib.Path format.

        Returns:
            tuple[Path]: A tuple of file paths.
        """
        return tuple(d.file for d in self)

    @property
    def names(self) -> tuple[str, ...]:
        """Returns the representative names of the files.

        Returns:
            tuple[str]: A tuple of file names.
        """
        return tuple(d.name for d in self)

    @property
    def checksums(self) -> dict[int, _FileWrapperABC]:
        """Returns a dictionary mapping sample CRC checksums to itself.

        Returns:
            dict[int, FileWrapperABC]: A dictionary of checksums and the object itself.
        """
        return {d.checksum: d for d in self}

    def load(self) -> None:
        """Loads data for all items in the list."""
        for d in self:
            d.load()

    def unload(self) -> None:
        """Unloads data for all items in the list."""
        for d in self:
            d.unload()

    def append(self, item: _FileWrapperABC) -> None:
        """Adds a new item to the list after validation.

        Args:
            item (FileWrapperABC): The item to append.

        Raises:
            TypeError: If the item is not of the correct type.
            KeyError: If the item already exists.
        """
        self._before_append_validate(item)
        self._data.append(item)

    def extend(self, other: Self) -> None:
        """Extends the list by appending items from another DataListABC object.

        Args:
            other (DataListABC): Another DataListABC object to extend from.
        """
        for item in other:
            self.append(item)

    def pop(self, i: int = -1) -> _FileWrapperABC:
        """Removes and returns the item at the given index.

        Args:
            i (int, optional): The index of the item to remove. Defaults to -1.

        Returns:
            FileWrapperABC: The removed item.
        """
        item: _FileWrapperABC = self._data.pop(i)
        return item

    def sort(self, /, *args, **kwargs) -> None:
        """Sorts the items in the list.

        Args:
            *args: Positional arguments for sorting.
            **kwargs: Keyword arguments for sorting.
        """
        self._data.sort(*args, **kwargs)

    def _before_append_validate(self, item: _FileWrapperABC) -> None:
        """Validates an item before appending it to the list.

        Args:
            item (FileWrapperABC): The item to validate.

        Raises:
            TypeError: If the item is not of the correct type.
            KeyError: If the item already exists in the list.
        """
        if not isinstance(item, self._bound_class):
            raise TypeError(
                f"{item.name}: Cannot add a {type(item).__qualname__} to a "
                f"{type(self).__qualname__} of {self._bound_class.__qualname__}."
            )
        if item in self:
            raise KeyError(f"{item.name}: Data already exists.")


class SeqDataListABC(DataListABC[_SeqFileWrapperABC]):
    """A list-like abstract base class providing partial list methods for managing SeqFileWrapperABC objects and their associated
    metadata.

    Attributes:
        files (tuple): A tuple of file paths.
        names (tuple): A tuple of file names.
        checksums (dict): A dictionary of names and their checksums.
    """

    __slots__ = ("_seqtype",)
    _bound_class: type[_SeqFileWrapperABC]

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(
        self, data: Sequence[str | Path | _SeqFileWrapperABC], *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO"
    ) -> None: ...
    @overload
    def __init__(
        self,
        data: Sequence[str | Path | _SeqFileWrapperABC],
        names: Sequence[str],
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    ) -> None: ...
    def __init__(
        self,
        data: Sequence[str | Path | _SeqFileWrapperABC] = (),
        names: Sequence[str] = (),
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    ) -> None:
        """Initializes the object and stores data into a list.

        Args:
            data (Sequence[str | Path | SeqFileWrapperABC], optional): A sequence of data items.
            names (Sequence[str], optional): A sequence of names corresponding to the data items.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Raises:
            RuntimeError: If names are provided but data is not.
            TypeError: If a data item cannot be converted to the bound class.
            KeyError: If the item already exists.
            SeqtypeError: If items represent different sequence types.
        """
        self._seqtype: Literal["dna", "pep", "rna", "NaN"] = "NaN"
        super().__init__(data, names, seqtype=seqtype)
        self._data: list[_SeqFileWrapperABC]

    def __repr__(self) -> str:
        """Returns a string representation of the object.

        Returns:
            str: The string representation of the object.
        """
        return type(self).__qualname__ + f"(seqtype={self.seqtype})" + "\n" + list_repr_wrapper(self._data)

    def __eq__(self, other: object) -> bool:
        """Checks if two objects have the same set of names.

        Args:
            other (SeqDataListABC): Another object to compare.

        Returns:
            bool: True if the sets of names are equal, False otherwise.

        Raises:
            SeqtypeError: If `self` and `other` represent different sequence types.
        """
        if not isinstance(other, type(self)):
            return NotImplemented
        check_seqtype(self, other)
        return super().__eq__(other)

    @property
    def seqtype(self) -> Literal["dna", "pep", "rna", "NaN"]:
        """Retrieves the common sequence type shared by all files in this object.

        Returns:
            str: The common sequence type of the files.
        """
        return self._seqtype

    def pop(self, i: int = -1) -> _SeqFileWrapperABC:
        """Removes and returns the item at the given index.

        Args:
            i (int, optional): The index of the item to remove. Defaults to -1.

        Returns:
            SeqFileWrapperABC: The removed item.
        """
        item = super().pop(i)
        if not self._data:
            self._seqtype = "NaN"
        return item

    def _before_append_validate(self, item: _SeqFileWrapperABC) -> None:
        """Validates an item before appending it to the list.

        Args:
            item (FileWrapperABC): The item to validate.

        Raises:
            TypeError: If the item is not of the correct type.
            KeyError: If the item already exists in the list.
            SeqtypeError: If `self` and `item` represent different sequence types.
        """
        super()._before_append_validate(item)
        seqtype = _add_seqtypes(self, item)
        self._seqtype = seqtype
