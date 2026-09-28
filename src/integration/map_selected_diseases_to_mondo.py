from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

MONDO_FILE = ROOT / "data/raw/mondo/mondo.obo"
OUT_FILE = ROOT / "data/interim/selected_diseases_mondo_mapping.csv"
REPORT_FILE = ROOT / "reports/selected_diseases_mondo_mapping_report.md"

SELECTED_DISEASES = [
    "epilepsy",
    "anxiety disorder",
    "stroke disorder",
    "Parkinson disease",
    "Alzheimer disease",
    "migraine disorder",
    "bipolar disorder",
    "schizophrenia",
    "dementia",
    "major depressive disorder",
    "attention deficit-hyperactivity disorder",
    "multiple sclerosis",
]


def flush_term(term, terms):
    if term.get("id", "").startswith("MONDO:"):
        terms.append(term)


def parse_mondo(path):
    terms = []
    current = {}

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")

            if line == "[Term]":
                if current:
                    flush_term(current, terms)
                current = {}
                continue

            if not line:
                continue

            if line.startswith("id: "):
                current["id"] = line[4:].strip()

            elif line.startswith("name: "):
                current["name"] = line[6:].strip()

            elif line.startswith("synonym: "):
                current.setdefault("synonyms", []).append(line)

        if current:
            flush_term(current, terms)

    return terms


def normalize(text):
    text = text.lower().strip()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text


def main():
    print("=" * 70)
    print("REMEDY - SELECTED DISEASE → CANONICAL MONDO MAPPING")
    print("=" * 70)

    print("\nParsing MONDO OBO...")

    terms = parse_mondo(MONDO_FILE)

    print(f"MONDO terms parsed: {len(terms):,}")

    # Exact label matches first
    rows = []

    normalized_lookup = {}

    for term in terms:
        name = term.get("name", "")

        if not name:
            continue

        normalized_lookup.setdefault(
            normalize(name),
            []
        ).append(term)

    print("\nSearching selected disease names...\n")

    for disease in SELECTED_DISEASES:
        matches = normalized_lookup.get(
            normalize(disease),
            []
        )

        print(f"{disease}")
        print("-" * len(disease))

        if not matches:
            print("NO EXACT LABEL MATCH")

            # Show loose candidates
            target = normalize(disease)
            candidates = []

            for term in terms:
                name = term.get("name", "")
                n = normalize(name)

                if (
                    target in n
                    or n in target
                    or any(
                        word in n.split()
                        for word in target.split()
                        if len(word) >= 5
                    )
                ):
                    candidates.append(term)

            candidates = candidates[:10]

            for candidate in candidates:
                print(
                    f"  {candidate.get('id')} | "
                    f"{candidate.get('name')}"
                )

        else:
            for match in matches:
                print(
                    f"  {match.get('id')} | "
                    f"{match.get('name')}"
                )

                rows.append({
                    "selected_disease": disease,
                    "mondo_id": match.get("id"),
                    "mondo_name": match.get("name"),
                    "match_type": "exact_label",
                })

        print()

    result = pd.DataFrame(rows)

    OUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    result.to_csv(
        OUT_FILE,
        index=False
    )

    report = [
        "# Selected Disease → MONDO Mapping Report",
        "",
        f"MONDO source: `{MONDO_FILE}`",
        "",
        "The script searches the MONDO release for exact normalized",
        "label matches for the 12 ReMedy target diseases.",
        "",
        "No MONDO identifier is guessed from external knowledge.",
        "",
        f"Exact matches found: {len(result):,}",
        "",
        "## Exact matches",
        "",
    ]

    if len(result):
        report.append(
            result.to_markdown(index=False)
        )
    else:
        report.append("No exact matches found.")

    report.extend([
        "",
        "## Important",
        "",
        "Disease mappings must be reviewed before being used in the",
        "final ReMedy knowledge graph.",
        "Grouped PrimeKG identifiers are not copied into MONDO IDs.",
    ])

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    REPORT_FILE.write_text(
        "\n".join(report),
        encoding="utf-8"
    )

    print("=" * 70)
    print("MONDO MAPPING SEARCH COMPLETE")
    print("=" * 70)
    print(f"Exact matches found: {len(result):,}")
    print(f"Saved: {OUT_FILE}")
    print(f"Saved: {REPORT_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()