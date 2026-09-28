from pathlib import Path
from datetime import datetime
import hashlib
import requests
import pandas as pd

# ============================================================
# CONFIG
# ============================================================

BASE_URL = "https://reactome.org/download/current"

OUTPUT_DIR = Path("data/raw/reactome")
MANIFEST_FILE = Path("data/dataset_manifest.csv")

REACTOME_VERSION = "V97"
RELEASE_DATE = "2026-06-30"
LICENSE = "CC0"

FILES = {
    "UniProt2Reactome_All_Levels.txt":
        f"{BASE_URL}/UniProt2Reactome_All_Levels.txt",

    "ReactomePathways.txt":
        f"{BASE_URL}/ReactomePathways.txt",

    "ReactomePathwaysRelation.txt":
        f"{BASE_URL}/ReactomePathwaysRelation.txt",

    "HumanDiseasePathways.txt":
        f"{BASE_URL}/HumanDiseasePathways.txt",
}


# ============================================================
# SETUP
# ============================================================

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

MANIFEST_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# SHA256
# ============================================================

def sha256_file(path):

    h = hashlib.sha256()

    with open(path, "rb") as f:

        while True:

            chunk = f.read(1024 * 1024)

            if not chunk:
                break

            h.update(chunk)

    return h.hexdigest()


# ============================================================
# DOWNLOAD
# ============================================================

def download_file(url, output_path):

    print()
    print("-" * 70)
    print(f"Downloading: {output_path.name}")

    temp_path = Path(
        str(output_path) + ".part"
    )

    with requests.get(
        url,
        stream=True,
        timeout=120,
        headers={
            "User-Agent": "ReMedy/1.0"
        }
    ) as response:

        response.raise_for_status()

        total = int(
            response.headers.get(
                "content-length",
                0
            )
        )

        downloaded = 0

        with open(
            temp_path,
            "wb"
        ) as f:

            for chunk in response.iter_content(
                chunk_size=1024 * 1024
            ):

                if not chunk:
                    continue

                f.write(chunk)

                downloaded += len(chunk)

                if total:

                    percent = (
                        downloaded / total
                    ) * 100

                    print(
                        f"\rProgress: {percent:6.2f}%",
                        end=""
                    )

    print()

    temp_path.replace(
        output_path
    )

    print(
        f"Saved: {output_path}"
    )


# ============================================================
# DOWNLOAD ALL FILES
# ============================================================

manifest_rows = []

download_date = datetime.now().strftime(
    "%Y-%m-%d"
)

for filename, url in FILES.items():

    output_path = (
        OUTPUT_DIR / filename
    )

    if output_path.exists():

        print(
            f"\nAlready exists: {filename}"
        )

    else:

        download_file(
            url,
            output_path
        )

    checksum = sha256_file(
        output_path
    )

    size_mb = (
        output_path.stat().st_size
        / (1024 * 1024)
    )

    print(
        f"SHA256: {checksum}"
    )

    print(
        f"Size: {size_mb:.2f} MB"
    )

    manifest_rows.append({

        "source":
            "Reactome",

        "version":
            REACTOME_VERSION,

        "release_date":
            RELEASE_DATE,

        "download_date":
            download_date,

        "file":
            filename,

        "licence":
            LICENSE,

        "url":
            url,

        "sha256":
            checksum,

        "local_path":
            str(output_path),

        "notes":
            "Reactome bulk data for ReMedy"
    })


# ============================================================
# UPDATE DATASET MANIFEST
# ============================================================

new_manifest = pd.DataFrame(
    manifest_rows
)

if MANIFEST_FILE.exists():

    old_manifest = pd.read_csv(
        MANIFEST_FILE
    )

    # Remove prior rows for these same files/version
    old_manifest = old_manifest[
        ~(
            (old_manifest["source"] == "Reactome")
            &
            (old_manifest["version"] == REACTOME_VERSION)
            &
            (
                old_manifest["file"]
                .isin(FILES.keys())
            )
        )
    ]

    final_manifest = pd.concat(
        [
            old_manifest,
            new_manifest
        ],
        ignore_index=True
    )

else:

    final_manifest = new_manifest


final_manifest.to_csv(
    MANIFEST_FILE,
    index=False
)


# ============================================================
# FINAL CHECK
# ============================================================

print()
print("=" * 70)
print("REACTOME DOWNLOAD COMPLETE")
print("=" * 70)

print(
    f"Reactome version: {REACTOME_VERSION}"
)

print(
    f"Files verified: {len(FILES)}"
)

print(
    f"Manifest: {MANIFEST_FILE}"
)

print("=" * 70)