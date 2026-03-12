"""Download module library."""

from __future__ import annotations

import csv
import hashlib
import shutil
import tarfile
import tempfile
from contextlib import ContextDecorator
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from . import BUSCO_URL, METADATA_FILE, logger


class BuscoParser(ContextDecorator):
    """Store/update local copy of HMM files into destination."""

    def __init__(self, *cfg_dirs: str | Path) -> None:
        """Initiate the object and download the latest metadata from online database."""
        self._cfg_dirs = cfg_dirs
        self._online_metadata: dict[str, dict] = {}
        self._local_metadata: dict[str, dict] = {}
        self._get_metadata_online()
        self._get_local_metadata()
        user_metadata = self._get_user_metadata()
        self._user_metadata_hash = hash(frozenset([(x, y["md5"]) for x, y in user_metadata.items()]))
        logger.debug(f"Original metadata hash: {self._user_metadata_hash}")

    def __enter__(self) -> BuscoParser:
        """Define the actions that will run when the object is created with `with` statement."""
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Define the actions that will run when the object is going to be destroyed by the end of the `with` statement."""
        logger.debug("Exiting context")
        user_metadata = self._get_user_metadata()
        user_metadata_hash = hash(frozenset([(x, y["md5"]) for x, y in user_metadata.items()]))
        logger.debug(f"Current metadata hash: {user_metadata_hash}")
        if user_metadata_hash != self._user_metadata_hash:
            logger.debug(f"Write changes to {self._cfg_dirs[0] / METADATA_FILE}")
            with open(self._cfg_dirs[0] / METADATA_FILE, "w") as f:
                writer = csv.writer(f, delimiter="\t")
                for name, info in user_metadata.items():
                    writer.writerow([name, info["md5"]])
        if exc_type:
            return False
        return True

    @property
    def online(self) -> list[str]:
        """Return a list of the markersets retrieve from online metadata."""
        return sorted(list(self._online_metadata))

    @property
    def local(self) -> list[str]:
        """Return a list of the markersets retrieve from local metadata."""
        local_markersets = []
        for markerset, info in self._local_metadata.items():
            if info["md5"] != self._online_metadata[markerset]["md5"]:
                markerset += " [Outdated]"
            local_markersets.append(markerset)
        return sorted(local_markersets)

    def close(self) -> None:
        self.__exit__(None, None, None)

    def download(self, markerset: str) -> None:
        """Download the given markerset from busco database."""
        if markerset in self._local_metadata:
            if self._local_metadata[markerset]["md5"] == self._online_metadata[markerset]["md5"]:
                logger.info("Markerset already exists and up to date.")
                return
            logger.info("Local markerset outdated. Update from online database...")
        else:
            logger.debug("Markerset not found. Download from online database...")
        data = _fetch_url(self._online_metadata[markerset]["url"])
        md5 = hashlib.md5(data).hexdigest()
        if md5 != self._online_metadata[markerset]["md5"]:
            raise HTTPError("The checksum of the downloaded contents do not match to the metadata.")
        self._save_markerset(data, markerset=markerset)
        self._local_metadata[markerset] = {"path": self._cfg_dirs[0], "md5": md5}

    def _get_metadata_online(self):
        """Get the metadata from busco url."""
        data = _fetch_url(f"{BUSCO_URL}/file_versions.tsv")
        for line in data.decode().split("\n"):
            line = line.split("\t")
            if line[-1] == "lineages":
                self._online_metadata[line[0]] = {
                    "url": f"{BUSCO_URL}/lineages/{line[0]}.{line[1]}.tar.gz",
                    "md5": line[2],
                }

    def _get_local_metadata(self):
        """Get the metadata from local file."""
        for cfg_dir in self._cfg_dirs[::-1]:  # Ensure the local metadata overwrite the global when overlapped
            metadata_file = cfg_dir / METADATA_FILE
            if metadata_file.is_file():
                with open(metadata_file) as f:
                    for line in csv.reader(f, delimiter="\t"):
                        if line[0].startswith("#"):
                            continue
                        if not (cfg_dir / line[0]).is_dir():
                            continue
                        self._local_metadata[line[0]] = {"path": cfg_dir, "md5": line[1]}
            else:
                logger.debug("Local metadata not found. Please download from the online database.")

    def _get_user_metadata(self):
        return {k: v for k, v in self._local_metadata.items() if v["path"] == self._cfg_dirs[0]} if self._local_metadata else {}

    def _save_markerset(self, data: bytes, markerset: str) -> None:
        """Save the content of the markerset to a local folder."""
        output = self._cfg_dirs[0] / markerset
        if output.is_dir():
            shutil.rmtree(output)
        output.mkdir(parents=True)
        with tempfile.NamedTemporaryFile("wb", delete=False) as fp:
            fp.write(data)
            fp.close()
            logger.debug("Extract files to %s.", {output})
            with tarfile.open(fp.name, "r:gz") as f:
                f.extractall(output.parent, filter="data")
        Path(fp.name).unlink()


def _fetch_url(url: str) -> bytes:
    """
    Fetch URL data content.

    Args:
        url (`str`): The url source.

    Returns:
        `bytes`: The content fetched from the url.
    """
    try:
        logger.debug("Download from %s ...", url)
        with urlopen(url) as response:
            content = response.read()
    except HTTPError as e:
        raise HTTPError("URL currently unavailable.") from e
    except URLError as e:
        raise URLError("URL not found or no internet connection available.") from e
    return content
