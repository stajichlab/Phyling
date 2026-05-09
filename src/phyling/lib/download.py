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
        self._cfg_dirs = [Path(x) for x in cfg_dirs]
        self._online_metadata: dict[str, dict] = {}
        self._local_metadata: dict[str, dict] = {}
        self._get_metadata_online()
        self._get_local_metadata()
        self._initial_hash = self._generate_metadata_hash()
        logger.debug("Original metadata hash %s", self._initial_hash)

    def __enter__(self) -> BuscoParser:
        """Define the actions that will run when the object is created with `with` statement."""
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        """Define the actions that will run when the object is going to be destroyed by the end of the `with` statement."""

        current_hash = self._generate_metadata_hash()
        logger.debug(f"Current metadata hash: {current_hash}")
        if current_hash != self._initial_hash:
            self._write_local_metadata()

        # Returning False allows exceptions to propagate normally
        return False

    @property
    def online(self) -> list[str]:
        """Return a list of the markersets retrieve from online metadata."""
        return sorted(list(self._online_metadata))

    @property
    def local(self) -> list[str]:
        """Return a list of the markersets retrieve from local metadata."""
        local_markersets = []
        for markerset, info in self._local_metadata.items():
            display_name = markerset
            if self._online_metadata and markerset in self._online_metadata:
                if info["md5"] != self._online_metadata[markerset]["md5"]:
                    display_name += " [Outdated]"
            local_markersets.append(display_name)
        return sorted(local_markersets)

    def close(self) -> None:
        self.__exit__(None, None, None)

    def download(self, markerset: str) -> None:
        """Download the given markerset from busco database."""
        if markerset not in self._online_metadata:
            raise KeyError(f"Markerset '{markerset}' not found in BUSCO database.")

        online_info = self._online_metadata[markerset]
        local_info = self._local_metadata.get(markerset)

        if local_info:
            if local_info["md5"] == online_info["md5"]:
                logger.info("Markerset %s already exists and is up to date.", markerset)
                return
            # Explicitly log the outdated status
            logger.info("Local markerset %s is outdated. Updating from online database...", markerset)
        else:
            logger.debug("Markerset %s not found locally. Downloading...", markerset)

        data = _fetch_url(online_info["url"])

        md5 = hashlib.md5(data).hexdigest()
        if md5 != online_info["md5"]:
            raise ValueError(f"MD5 mismatch for {markerset}. Expected {online_info['md5']}, got {md5}")

        self._save_markerset(data, markerset)
        self._local_metadata[markerset] = {"path": self._cfg_dirs[0], "md5": md5}

    def _get_metadata_online(self):
        """Get the metadata from busco url."""
        try:
            data = _fetch_url(f"{BUSCO_URL}/file_versions.tsv")
            for line in data.decode().splitlines():
                parts = line.split("\t")
                if parts[-1] == "lineages":
                    self._online_metadata[parts[0]] = {
                        "url": f"{BUSCO_URL}/lineages/{parts[0]}.{parts[1]}.tar.gz",
                        "md5": parts[2],
                    }
        except Exception as e:
            logger.error("Could not retrieve online metadata: %s", e)

    def _get_local_metadata(self):
        """Get the metadata from local file."""
        for cfg_dir in self._cfg_dirs[::-1]:  # Ensure the local metadata overwrite the global when overlapped
            metadata_file = cfg_dir / METADATA_FILE
            if not metadata_file.is_file():
                continue
            with open(metadata_file) as f:
                for line in csv.reader(f, delimiter="\t"):
                    if line[0].startswith("#"):
                        continue
                    if (cfg_dir / line[0]).is_dir():
                        self._local_metadata[line[0]] = {"path": cfg_dir, "md5": line[1]}

    def _write_local_metadata(self) -> None:
        """Atomically write metadata to the primary config directory."""
        dest = self._cfg_dirs[0] / METADATA_FILE
        dest.parent.mkdir(parents=True, exist_ok=True)

        logger.debug("Writing changes to %s", dest)
        # Write to temp file first then rename to prevent corruption on crash
        with tempfile.NamedTemporaryFile("w", dir=dest.parent, delete=False) as f:
            writer = csv.writer(f, delimiter="\t")
            for name, info in self._get_user_metadata().items():
                writer.writerow([name, info["md5"]])
        Path(f.name).replace(dest)

    def _get_user_metadata(self):
        return {k: v for k, v in self._local_metadata.items() if v["path"] == self._cfg_dirs[0]}

    def _generate_metadata_hash(self) -> int:
        """Generate a hash of the current local metadata for change detection."""
        user_data = self._get_user_metadata()
        # Sort keys to ensure deterministic hashing
        state = tuple(sorted((k, v["md5"]) for k, v in user_data.items()))
        return hash(state)

    def _save_markerset(self, data: bytes, markerset: str) -> None:
        """Save the content of the markerset to a local folder."""
        output_dir = self._cfg_dirs[0] / markerset

        # 1. Use a temporary directory for safe extraction
        with tempfile.TemporaryDirectory(dir=self._cfg_dirs[0]) as tmp_dir:
            tmp_path = Path(tmp_dir)
            tar_path = tmp_path / "data.tar.gz"
            tar_path.write_bytes(data)

            with tarfile.open(tar_path, "r:gz") as f:
                # Python 3.12+ Security: 'data' filter prevents path traversal
                if hasattr(tarfile, "data_filter"):
                    f.extractall(tmp_path, filter="data")
                else:
                    # Python 3.9 - 3.11 behavior
                    f.extractall(tmp_path)

            # 2. Cleanup existing and move new content into place
            if output_dir.exists():
                shutil.rmtree(output_dir)

            # 3. Find the extracted folder
            # BUSCO tarballs usually contain a single top-level folder
            extracted_folders = [p for p in tmp_path.iterdir() if p.is_dir()]
            if extracted_folders:
                # Move the first directory found to the final destination
                shutil.move(str(extracted_folders[0]), str(output_dir))
            else:
                raise RuntimeError(f"Failed to find extracted folder for {markerset}")


def _fetch_url(url: str, timeout: int = 30) -> bytes:
    """
    Fetch URL data content.

    Args:
        url (`str`): The url source.

    Returns:
        `bytes`: The content fetched from the url.
    """
    try:
        logger.debug("Download from %s ...", url)
        with urlopen(url, timeout=timeout) as response:
            return response.read()
    except HTTPError as e:
        raise HTTPError(url, e.code, f"Server returned error: {e.reason}", e.hdrs, None) from e
    except URLError as e:
        raise URLError(f"Connection failed to {url}: {e.reason}") from e
