from __future__ import annotations

import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, TextIO

import numpy as np

from ..lib import TreeMethods
from ..lib._abc import list_repr_wrapper

# Model name in FastTree, Raxml-ng, IQ-Tree
DNA_MODELS = np.array(
    [
        ("GTR", "GTR", "GTR"),
        ("", "SYM", "SYM"),
        ("", "TVM", "TVM"),
        ("", "TVMef", "TVMe"),
        ("", "TIM3uf", "TIM3"),
        ("", "TIM3", "TIM3e"),
        ("", "TIM2uf", "TIM2"),
        ("", "TIM2", "TIM2e"),
        ("", "TIM1uf", "TIM"),
        ("", "TIM1", "TIMe"),
        ("", "TPM3uf", "TPM3u"),
        ("", "TPM3", "TPM3"),
        ("", "TPM2uf", "TPM2u"),
        ("", "TPM2", "TPM2"),
        ("", "K81uf", "K3Pu"),
        ("", "K81", "K3P"),
        ("", "TN93", "TN"),
        ("", "TN93ef", "TNe"),
        ("", "HKY", "HKY"),
        ("", "K80", "K2P"),
        ("", "F81", "F81"),
        ("JC", "JC", "JC"),
    ]
)


PEP_MODELS = np.array(
    [
        ("LG", "LG", "LG"),
        ("WAG", "WAG", "WAG"),
        ("JTT", "JTT", "JTT"),
        ("", "Q.PFAM", "Q.PFAM"),
        ("", "Q.BIRD", "Q.BIRD"),
        ("", "Q.MAMMAL", "Q.MAMMAL"),
        ("", "Q.INSECT", "Q.INSECT"),
        ("", "Q.PLANT", "Q.PLANT"),
        ("", "Q.YEAST", "Q.YEAST"),
        ("", "JTT-DCMUT", "JTTDCMUT"),
        ("", "DCMUT", "DCMUT"),
        ("", "VT", "VT"),
        ("", "PMB", "PMB"),
        ("", "BLOSUM62", "BLOSUM62"),
        ("", "DAYHOFF", "DAYHOFF"),
        ("", "MTREV", "MTREV"),
        ("", "MTART", "MTART"),
        ("", "MTZOA", "MTZOA"),
        ("", "", "MTMET"),
        ("", "", "mtVer"),
        ("", "", "MTINV"),
        ("", "MTMAM", "MTMAM"),
        ("", "", "FLAVI"),
        ("", "HIVB", "HIVB"),
        ("", "HIVW", "HIVW"),
        ("", "FLU", "FLU"),
        ("", "RTREV", "RTREV"),
        ("", "CPREV", "CPREV"),
    ]
)

ALL_MODELS = np.concat([DNA_MODELS, PEP_MODELS])

STATIONARY_CODES = np.array(
    [
        ("", "F", "F"),
        ("", "FC", "F"),
        ("", "FU", "F"),
        ("", "FE", "FQ"),
        ("", "FO", "FO"),
    ]
)

INVARIANT_CODES = np.array(
    [
        ("", "I", "I"),
        ("", "IO", "I"),
        ("", "IC", "I"),
        ("", "IU", "I"),
    ]
)


@dataclass
class PartitionRecord:
    """A generalized container for partitioning info used in RAxML-NG and IQTree partitioning analyses."""

    name: str
    start: int
    end: int
    model: str = ""
    stationary: str = ""
    invariant: str = ""
    gamma: str = ""


class Partitions(ABC):
    """A collection of PartitionRecord objects used in RAxML-NG and IQTree partitioning analyses.

    This class acts as a container for PartitionRecord instances, allowing seamless conversion between formats used by RAxML-NG
    and IQTree.
    """

    _model_tbl: np.ndarray = ALL_MODELS
    _stationary_tbl: np.ndarray = STATIONARY_CODES
    _invariant_tbl: np.ndarray = INVARIANT_CODES
    _methods: tuple[str] = tuple(x.name.lower() for x in TreeMethods)

    __slots__ = ("_data", "_method_idx", "_avail_models")

    def __init__(self, source: Literal["ft", "raxml", "iqtree"]):
        """Initialize a Partitions object.

        Args:
            source (Literal["ft", "raxml", "iqtree"]): The format of this Partitions object.
                - "ft": FastTree
                - "raxml": RAxML-NG
                - "iqtree": IQTree
        """
        self._data: list[PartitionRecord] = []
        self._method_idx: int = self._methods.index(source)
        self._avail_models: np.ndarray = self._model_tbl[self._model_tbl[:, self._method_idx] != "", self._method_idx]

    def __repr__(self):
        """Return a string representation of the Partitions.

        Returns:
            str: The string representation of the Partitions.
        """
        return type(self).__qualname__ + f"(seqtype={self.seqtype})" + "\n" + list_repr_wrapper(self._data)

    def __iter__(self):
        """Iterate over the PartitionRecord.

        Returns:
            Iterator[PartitionRecord]: An iterator over the PartitionRecord.
        """
        return iter(self._data)

    def __len__(self) -> int:
        """Get the number of PartitionRecord.

        Returns:
            int: The number of PartitionRecord.
        """
        return len(self._data)

    def __getitem__(self, key: int) -> PartitionRecord:
        """Retrieves an PartitionRecord by index.

        Args:
            key (int): The index to retrieve.

        Returns:
            PartitionRecord: The corresponding PartitionRecord.
        """
        if isinstance(key, int):
            return self._data[key]

    def add(self, part_rec: PartitionRecord):
        """Add a new search PartitionRecord to the Partitions.

        Args:
            part_rec (PartitionRecord): The PartitionRecord to add.
        """
        model, _ = self._decipher_param(part_rec.model)
        if model not in self._avail_models:
            raise RuntimeError(f"Model {model} is not supported. Supported models: {', '.join(self._avail_models)}")
        if part_rec.name in (msa.name for msa in self):
            raise RuntimeError(f"Partition name '{part_rec.name}' already exists.")
        self._data.append(part_rec)

    def convert_to(self, format: Literal["raxml", "iqtree"]) -> Partitions:
        """Convert the model, stationary, and invariant codes of each PartitionRecord in the Partitions to the specified format.

        Args:
            format (Literal["raxml", "iqtree"]): The target format.
                - "raxml": RAxML-NG
                - "iqtree": IQTree

        Returns:
            Partitions: A new Partitions object with codes converted to the specified format.
        """

        def decipher_gamma(gamma):
            return re.sub(r"^(G\d+)m", lambda m: m.group(1), gamma)

        if format not in ("iqtree", "raxml"):
            raise RuntimeError(f"Partitioning is not supported in {format}.")

        new_parts = type(self)(source=format)
        to_method_idx = self._methods.index(format)
        for msa in self:
            model, p = self._decipher_param(msa.model)
            model = self._model_tbl[self._model_tbl[:, self._method_idx] == model, to_method_idx].item()
            if not model:
                raise RuntimeError(f"Model {msa.model} is not mappable to {format} format.")
            if p:
                model += f"{{{','.join([str(p_) for p_ in p])}}}"
            if msa.stationary:
                stationary, p = self._decipher_param(msa.stationary)
                stationary = self._stationary_tbl[self._stationary_tbl[:, self._method_idx] == stationary, to_method_idx][
                    0
                ].item()
                if not stationary:
                    raise RuntimeError(f"Stationary code {msa.stationary} is not mappable to {format} format.")
                if p:
                    stationary += f"{{{','.join([str(p_) for p_ in p])}}}"
            if msa.invariant:
                invariant, p = self._decipher_param(msa.invariant)
                invariant = self._invariant_tbl[self._invariant_tbl[:, self._method_idx] == invariant, to_method_idx][0].item()
                if not invariant:
                    raise RuntimeError(f"Invariant code {msa.invariant} is not mappable to {format} format.")
                if p:
                    invariant += f"{{{','.join([str(p_) for p_ in p])}}}"
            if msa.gamma:
                gamma = decipher_gamma(msa.gamma)

            new_parts.add(
                PartitionRecord(
                    name=msa.name,
                    start=msa.start,
                    end=msa.end,
                    model=model,
                    stationary=stationary if msa.stationary else "",
                    invariant=invariant if msa.invariant else "",
                    gamma=gamma if msa.gamma else "",
                )
            )
        return new_parts

    def _decipher_param(self, param: str):
        name, values = re.match(r"^([^{]+)(?:\{([^}]+)\})?$", param).groups()
        if values:
            values = [float(v) for v in values.split(",")]
        return name, values


class CustomHandler(ABC):
    """A custom text file handler abstract class."""

    __slots__ = ("_file", "_format")

    def __init__(self, file: str | Path, mode: str = "r", *, format: Literal["raxml", "iqtree"]):
        """Initialize a CustomHandler object.

        Args:
            file (str | Path): The path of the file.
            mode (str): 'r' for reading and 'w' for writing.
            format (Literal["raxml", "iqtree"]): The format (model, stationary and invariant codes) uses in the file.
                - "raxml": RAxML-NG
                - "iqtree": IQTree
        """
        if "b" in mode:
            raise ValueError("File must be opened in text mode, not binary mode.")
        self._file: TextIO = open(file, mode)
        self._format = format

    def __enter__(self):
        """Define the actions that will run when the object is created with `with` statement."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Define the actions that will run when the object is going to be destroyed by the end of the `with` statement."""
        self.close()

    def close(self):
        """Closes the file."""
        self._file.close()

    @abstractmethod
    def read(self): ...

    @abstractmethod
    def write(self): ...

    def readable(self):
        return self._file.readable()

    def writable(self):
        return self._file.writable()

    def seek(self, offset: int, whence: int = 0, /):
        self._file.seek(offset, whence)


class NexusHandler(CustomHandler):
    """A subclass of CustomHandler to help parse the nexus partition file."""

    _header: str = "#nexus"
    _delimiter: str = ";"

    def __init__(self, file: str | Path, mode: str = "r", *, format: Literal["raxml", "iqtree"] = "iqtree"):
        """Initialize a NexusHandler object.

        Args:
            file (str | Path): The path of the file.
            mode (str): 'r' for reading and 'w' for writing.
            format (Literal["raxml", "iqtree"]): The format (model, stationary and invariant codes) uses in the file.
                - "raxml": RAxML-NG
                - "iqtree": IQTree
        """
        super().__init__(file, mode, format=format)
        if "r" in mode and self._file.readline().strip().lower() != self._header:
            raise RuntimeError("Not a nexus file.")

    def read(self) -> dict[str, Partitions]:
        """Reads from the file.

        Currently only supports "sets" block reading.
        """
        if self._file.tell() < 7:
            self._file.seek(7)
        r = {}
        blocks = self._parse_block()
        while blocks:
            try:
                title, blocklines = next(blocks)
            except StopIteration:
                break
            try:
                r[title] = getattr(self, f"_parse_{title}")(blocklines)
            except AttributeError:
                raise RuntimeWarning(f"{title} block processing is not implemented yet.")
        return r

    def write(self, data: dict[str, Partitions]):
        """Writes data to the file."""
        self._file.write(f"{self._header}\n")
        for block in data:
            self._file.write(f"begin {block}{self._delimiter}\n")
            charpartition = []
            for partition in data[block]:
                self._file.write(f"  charset {partition.name} = {partition.start}-{partition.end}{self._delimiter}\n")
                part_info = [partition.model]
                for attr in ("stationary", "invariant", "gamma"):
                    v = getattr(partition, attr)
                    if v:
                        part_info.append(v)
                charpartition.append(f"    {'+'.join(part_info)}: {partition.name}")
            charpartition = "  charpartition mymodel =\n" + ",\n".join(charpartition)
            self._file.write(f"{charpartition};\n")
            self._file.write("end;\n")

    def _parse_line(self):
        """Reads up to the first occurrence of the custom delimiter."""
        line = ""
        while True:
            char = self._file.read(1)
            if not char or char == self._delimiter:  # Stop at delimiter or EOF
                break
            line += char
        return line.strip()

    def _parse_block(self):
        inblock = False
        blocklines = []
        line = self._parse_line()
        while line:
            if line.lower().startswith("begin"):
                if not inblock:
                    inblock = True
                    title = line.split()[1].lower()
                else:
                    raise RuntimeError(f"Illegal block nesting in block {title}")
            elif line.lower().startswith("end"):
                if inblock:
                    inblock = False
                    yield title, blocklines
                    blocklines.clear()
                else:
                    raise RuntimeError("Unmatched 'end'.")
            elif inblock:
                blocklines.append(line)
            line = self._parse_line()

    def _parse_sets(self, blocklines: list[str]):
        partitions: dict[str, PartitionRecord] = {}
        while blocklines:
            line = blocklines.pop(0)
            line = line.lstrip("\n")
            if line.startswith("charset"):
                part_name, content = line[7:].replace(" ", "").split("=")
                start, end = [int(i) for i in content.split(",")[-1].split("-")]
                partitions[part_name] = PartitionRecord(name=part_name, start=start, end=end)
            if line.startswith("charpartition"):
                _, content = line[13:].replace(" ", "").split("=")
                for charpartition in content.strip().split(",\n"):
                    params, part_name = charpartition.split(":")
                    part_name = re.sub(r"\{[^}]*\}$", "", part_name)
                    params = params.split("+")
                    partitions[part_name].model = params.pop(0)
                    while params:
                        param = params.pop(0)
                        if param.startswith("F"):
                            partitions[part_name].stationary = param
                        elif param.startswith("I"):
                            partitions[part_name].invariant = param
                        elif param.startswith("G"):
                            partitions[part_name].gamma = param
        partition_obj = Partitions(self._format)

        # return partitions
        for p in partitions.values():
            partition_obj.add(p)
        return partition_obj


class RaxmlHandler(CustomHandler):
    """A subclass of CustomHandler to help parse the raxml-style partition file."""

    _delimiter: str = "\n"

    def __init__(self, file: str | Path, mode: str = "r", *, format: Literal["raxml", "iqtree"] = "raxml"):
        """Initialize a RaxmlHandler object.

        Args:
            file (str | Path): The path of the file.
            mode (str): 'r' for reading and 'w' for writing.
            format (Literal["raxml", "iqtree"]): The format (model, stationary and invariant codes) uses in the file.
                - "raxml": RAxML-NG
                - "iqtree": IQTree
        """
        super().__init__(file, mode, format=format)

    def read(self):
        """Reads from the file."""
        partition_obj = Partitions(self._format)
        while True:
            line = self._parse_line()
            if not line:
                break
            params, part_info = line.replace(" ", "").split(",")
            part_name, content = part_info.split("=")
            start, end = [int(i) for i in content.split("-")]
            p = PartitionRecord(name=part_name, start=start, end=end)
            params = re.sub(r"\+BU\{[^}]+\}", "", params)
            params = re.sub(r"(\{[^}]+\})", lambda m: m.group(1).replace("/", ","), params)
            params = params.split("+")
            p.model = params.pop(0)
            while params:
                param = params.pop(0)
                if param.startswith("F"):
                    p.stationary = param
                elif param.startswith("I"):
                    p.invariant = param
                elif param.startswith("G"):
                    p.gamma = param
            partition_obj.add(p)
        return partition_obj

    def write(self, data: Partitions, param_sep="/"):
        """Writes data to the file."""
        for partition in data:
            params = [partition.model]
            for attr in ("stationary", "invariant", "gamma"):
                v = getattr(partition, attr)
                if v:
                    params.append(v)
            params = re.sub(r"(\{[^}]+\})", lambda m: m.group(1).replace(",", param_sep), "+".join(params))
            self._file.write(f"{params}, {partition.name} = {partition.start}-{partition.end}{self._delimiter}")

    def _parse_line(self):
        """Reads up to the first occurrence of the custom delimiter."""
        line = ""
        while True:
            char = self._file.read(1)
            if not char or char == self._delimiter:  # Stop at delimiter or EOF
                break
            line += char
        return line.strip()
