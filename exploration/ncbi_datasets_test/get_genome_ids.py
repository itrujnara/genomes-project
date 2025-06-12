import sys
import json


def get_genome_ids(json_file: str) -> None:
    with open(json_file, "r") as f:
        data = json.load(f)
        for genome in data["reports"]:
            print(genome["accession"])


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python get_genome_ids.py <json>")
        sys.exit(1)
    get_genome_ids(sys.argv[1])


if __name__ == "__main__":
    main()