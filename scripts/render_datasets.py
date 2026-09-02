#!/usr/bin/env python3
"""Validate datasets.yml and render its table into README.md."""

from __future__ import annotations

import argparse
import difflib
import sys
from pathlib import Path, PurePosixPath
from urllib.parse import quote, urlparse

try:
    import yaml
except ImportError as exc:  # pragma: no cover - depends on the local environment
    raise SystemExit(
        "PyYAML is required. Install it with "
        "'python -m pip install -r scripts/requirements.txt'."
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "datasets.yml"
README = ROOT / "README.md"
START_MARKER = "<!-- DATASETS:START -->"
END_MARKER = "<!-- DATASETS:END -->"
REQUIRED_DATASET_FIELDS = {
    "id",
    "name",
    "data_type",
    "source_url",
    "analyses",
}
OPTIONAL_DATASET_FIELDS = {
    "version",
    "citation",
    "license",
    "license_url",
    "downloader",
}
ANALYSIS_FIELDS = {"label", "path"}


class ManifestError(ValueError):
    """Raised when datasets.yml does not satisfy the repository contract."""


def fail(message: str) -> None:
    raise ManifestError(message)


def require_text(value: object, location: str) -> str:
    if not isinstance(value, str) or not value.strip():
        fail(f"{location} must be a non-empty string")
    if "\n" in value or "|" in value:
        fail(f"{location} cannot contain a newline or table separator ('|')")
    return value.strip()


def validate_repo_path(value: object, location: str) -> str:
    path_text = require_text(value, location)
    posix_path = PurePosixPath(path_text)
    if posix_path.is_absolute() or ".." in posix_path.parts:
        fail(f"{location} must be a repository-relative path without '..'")

    resolved = (ROOT / Path(*posix_path.parts)).resolve()
    if ROOT not in resolved.parents or not resolved.is_file():
        fail(f"{location} does not point to an existing repository file: {path_text}")
    return path_text


def load_and_validate_manifest() -> list[dict[str, object]]:
    try:
        document = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        fail(f"datasets.yml is not valid YAML: {exc}")

    if not isinstance(document, dict):
        fail("datasets.yml must contain a mapping at its root")
    if set(document) != {"schema_version", "datasets"}:
        fail("datasets.yml must contain exactly 'schema_version' and 'datasets'")
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        fail("unsupported schema_version; expected 1")

    datasets = document["datasets"]
    if not isinstance(datasets, list) or not datasets:
        fail("datasets must be a non-empty list")

    seen_ids: set[str] = set()
    for index, dataset in enumerate(datasets, start=1):
        location = f"datasets[{index}]"
        if not isinstance(dataset, dict):
            fail(f"{location} must be a mapping")

        fields = set(dataset)
        missing = REQUIRED_DATASET_FIELDS - fields
        unknown = fields - REQUIRED_DATASET_FIELDS - OPTIONAL_DATASET_FIELDS
        if missing:
            fail(f"{location} is missing: {', '.join(sorted(missing))}")
        if unknown:
            fail(f"{location} has unknown fields: {', '.join(sorted(unknown))}")

        dataset_id = require_text(dataset["id"], f"{location}.id")
        if dataset_id.casefold() in seen_ids:
            fail(f"duplicate dataset id: {dataset_id}")
        seen_ids.add(dataset_id.casefold())

        require_text(dataset["name"], f"{location}.name")
        require_text(dataset["data_type"], f"{location}.data_type")

        source_url = require_text(dataset["source_url"], f"{location}.source_url")
        parsed_url = urlparse(source_url)
        if parsed_url.scheme != "https" or not parsed_url.netloc:
            fail(f"{location}.source_url must be an absolute HTTPS URL")

        for optional_text in ("version", "citation", "license"):
            if optional_text in dataset:
                require_text(dataset[optional_text], f"{location}.{optional_text}")
        if "license_url" in dataset:
            if "license" not in dataset:
                fail(f"{location}.license_url requires a license field")
            license_url = require_text(
                dataset["license_url"], f"{location}.license_url"
            )
            parsed_license_url = urlparse(license_url)
            if parsed_license_url.scheme != "https" or not parsed_license_url.netloc:
                fail(f"{location}.license_url must be an absolute HTTPS URL")
        if "downloader" in dataset:
            validate_repo_path(dataset["downloader"], f"{location}.downloader")

        analyses = dataset["analyses"]
        if not isinstance(analyses, list) or not analyses:
            fail(f"{location}.analyses must be a non-empty list")
        seen_analysis_paths: set[str] = set()
        for analysis_index, analysis in enumerate(analyses, start=1):
            analysis_location = f"{location}.analyses[{analysis_index}]"
            if not isinstance(analysis, dict) or set(analysis) != ANALYSIS_FIELDS:
                fail(f"{analysis_location} must contain exactly 'label' and 'path'")
            require_text(analysis["label"], f"{analysis_location}.label")
            analysis_path = validate_repo_path(
                analysis["path"], f"{analysis_location}.path"
            )
            if analysis_path.casefold() in seen_analysis_paths:
                fail(f"{location} contains a duplicate analysis path: {analysis_path}")
            seen_analysis_paths.add(analysis_path.casefold())

    return datasets


def markdown_path(path_text: str) -> str:
    return quote(path_text, safe="/._-")


def render_table(datasets: list[dict[str, object]]) -> str:
    lines = [
        START_MARKER,
        "<!-- Generated from datasets.yml by scripts/render_datasets.py. Do not edit manually. -->",
        "",
        "| Data type | Dataset | Used in |",
        "| --- | --- | --- |",
    ]
    for dataset in datasets:
        display_id = str(dataset["id"])
        if "version" in dataset:
            display_id += f" v{dataset['version']}"
        dataset_cell = f"[{display_id}]({dataset['source_url']}) — {dataset['name']}"
        if "citation" in dataset:
            dataset_cell += f" ({dataset['citation']})"
        if "license" in dataset:
            license_text = str(dataset["license"])
            if "license_url" in dataset:
                license_text = f"[{license_text}]({dataset['license_url']})"
            dataset_cell += f" — License: {license_text}"

        analysis_links = []
        for analysis in dataset["analyses"]:  # type: ignore[union-attr]
            analysis_links.append(
                f"[{analysis['label']}]({markdown_path(str(analysis['path']))})"
            )
        lines.append(
            f"| {dataset['data_type']} | {dataset_cell} | {', '.join(analysis_links)} |"
        )

    lines.extend([END_MARKER, ""])
    return "\n".join(lines)


def update_readme(table: str) -> tuple[str, str]:
    current = README.read_text(encoding="utf-8")
    if current.count(START_MARKER) != 1 or current.count(END_MARKER) != 1:
        fail("README.md must contain exactly one dataset start marker and one end marker")
    start = current.index(START_MARKER)
    end = current.index(END_MARKER, start) + len(END_MARKER)
    return current, current[:start] + table.rstrip("\n") + current[end:]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if README.md is not synchronized; do not write files",
    )
    args = parser.parse_args()

    try:
        datasets = load_and_validate_manifest()
        current, expected = update_readme(render_table(datasets))
    except (ManifestError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if current == expected:
        print("README dataset table is up to date.")
        return 0

    if args.check:
        print(
            "README dataset table is out of date. Run scripts/render_datasets.py.",
            file=sys.stderr,
        )
        sys.stderr.writelines(
            difflib.unified_diff(
                current.splitlines(keepends=True),
                expected.splitlines(keepends=True),
                fromfile="README.md",
                tofile="README.md (generated)",
            )
        )
        return 1

    README.write_text(expected, encoding="utf-8", newline="\n")
    print(f"Updated README.md with {len(datasets)} datasets.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
