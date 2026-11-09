import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


PARSER = Path("src/parse_pooled_inner_zip.py")

TMP_DIR = Path("tmp/pooled_runner")
PARTS_DIR = Path("results/pooled_parts")
STATE_DIR = Path("results/pooled_state")

COMPLETED_FILE = STATE_DIR / "completed_inner_zips.txt"
FAILED_FILE = STATE_DIR / "failed_inner_zips.tsv"


def safe_name(member_name):
    """
    Convert an inner ZIP path to a safe output filename.

    Example:
    AF3_completed_zips/folds_2025_03_18_20_18.zip
    ->
    folds_2025_03_18_20_18
    """

    name = Path(member_name).name

    if name.lower().endswith(".zip"):
        name = name[:-4]

    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        name,
    )


def load_completed():
    """Load successfully completed inner ZIP members."""

    if not COMPLETED_FILE.exists():
        return set()

    with open(COMPLETED_FILE) as handle:
        return {
            line.strip()
            for line in handle
            if line.strip()
        }


def mark_completed(member_name):
    """Record a successfully processed inner ZIP."""

    with open(
        COMPLETED_FILE,
        "a",
    ) as handle:
        handle.write(
            member_name + "\n"
        )


def record_failure(
    member_name,
    return_code,
):
    """Record a failed inner ZIP."""

    with open(
        FAILED_FILE,
        "a",
    ) as handle:
        handle.write(
            f"{member_name}\t{return_code}\n"
        )


def find_inner_zips(outer_zip):
    """Find real inner ZIP files, ignoring macOS metadata."""

    members = []

    for info in outer_zip.infolist():

        name = info.filename

        if info.is_dir():
            continue

        if not name.lower().endswith(".zip"):
            continue

        if "__MACOSX/" in name:
            continue

        if Path(name).name.startswith("._"):
            continue

        members.append(info)

    return members


def main(outer_zip_path):

    outer_zip_path = Path(
        outer_zip_path
    )

    if not outer_zip_path.exists():
        raise FileNotFoundError(
            outer_zip_path
        )

    TMP_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    PARTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    STATE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    completed = load_completed()

    print(
        f"Opening outer ZIP: "
        f"{outer_zip_path}"
    )

    with zipfile.ZipFile(
        outer_zip_path
    ) as outer:

        inner_zips = find_inner_zips(
            outer
        )

        print(
            f"Found {len(inner_zips)} "
            f"inner ZIP files."
        )

        for index, info in enumerate(
            inner_zips,
            start=1,
        ):

            member_name = info.filename
            stem = safe_name(
                member_name
            )

            output_csv = (
                PARTS_DIR
                / f"{stem}.csv"
            )

            # Only skip if both the state record
            # and the result file exist.
            if (
                member_name in completed
                and output_csv.exists()
            ):
                print(
                    f"[{index}/{len(inner_zips)}] "
                    f"SKIP completed: "
                    f"{member_name}"
                )

                continue

            print()
            print(
                "=" * 70
            )

            print(
                f"[{index}/{len(inner_zips)}] "
                f"{member_name}"
            )

            print(
                f"Compressed size: "
                f"{info.compress_size / 1024**3:.2f} GB"
            )

            print(
                f"Uncompressed inner ZIP size: "
                f"{info.file_size / 1024**3:.2f} GB"
            )

            temp_zip = (
                TMP_DIR
                / f"{stem}.zip"
            )

            # Remove stale temporary file
            # from a previous interrupted run.
            if temp_zip.exists():
                print(
                    f"Removing stale temp file: "
                    f"{temp_zip}"
                )

                temp_zip.unlink()

            # ---------------------------------
            # Extract ONE inner ZIP
            # ---------------------------------

            print(
                "Extracting inner ZIP..."
            )

            with outer.open(
                info
            ) as source:

                with open(
                    temp_zip,
                    "wb",
                ) as destination:

                    shutil.copyfileobj(
                        source,
                        destination,
                        length=16 * 1024 * 1024,
                    )

            print(
                f"Temporary file ready: "
                f"{temp_zip}"
            )

            # ---------------------------------
            # Run pooled parser
            # ---------------------------------

            print(
                "Parsing..."
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(PARSER),
                    str(temp_zip),
                    str(output_csv),
                ]
            )

            # ---------------------------------
            # Handle success / failure
            # ---------------------------------

            if result.returncode == 0:

                mark_completed(
                    member_name
                )

                completed.add(
                    member_name
                )

                print(
                    f"SUCCESS: {member_name}"
                )

                print(
                    f"Output: {output_csv}"
                )

            else:

                record_failure(
                    member_name,
                    result.returncode,
                )

                print(
                    f"FAILED: {member_name}",
                    file=sys.stderr,
                )

                print(
                    f"Parser exit code: "
                    f"{result.returncode}",
                    file=sys.stderr,
                )

                print(
                    "Stopping so the error can "
                    "be inspected.",
                    file=sys.stderr,
                )

                if temp_zip.exists():
                    temp_zip.unlink()

                sys.exit(
                    result.returncode
                )

            # ---------------------------------
            # Delete temporary inner ZIP
            # ---------------------------------

            if temp_zip.exists():
                temp_zip.unlink()

            print(
                "Temporary inner ZIP deleted."
            )

    print()
    print(
        "=" * 70
    )

    print(
        "All inner ZIP files completed."
    )


if __name__ == "__main__":

    if len(sys.argv) != 2:

        print(
            "Usage: python "
            "run_pooled_outer_zip.py "
            "<outer_zip>"
        )

        sys.exit(1)

    main(
        sys.argv[1]
    )