import sys
import json


def get_species(json_file: str) -> None:
    with open(json_file, "r") as f:
        data = json.load(f)
        for species in data["reports"]:
            print(species["taxonomy"]["tax_id"])


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python get_species.py <json>")
        sys.exit(1)
    get_species(sys.argv[1])


if __name__ == "__main__":
    main()