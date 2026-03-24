"""CLI entry point for stangene."""

import argparse
import sys

from stangene._logging import get_logger

logger = get_logger("cli")


def main():
    parser = argparse.ArgumentParser(
        prog="stangene",
        description="Gene identifier harmonization for single-cell transcriptomics",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    harm_parser = subparsers.add_parser("harmonize", help="Harmonize gene identifiers")
    harm_parser.add_argument("--input", required=True, help="Path to input file")
    harm_parser.add_argument("--species", required=True, help="Species name")
    harm_parser.add_argument("--output-dir", required=True, help="Output directory")
    harm_parser.add_argument("--dataset-name", default=None, help="Dataset name")
    harm_parser.add_argument("--reference-dir", default=None, help="Reference directory")

    build_parser = subparsers.add_parser("build-refs", help="Build reference databases")
    build_parser.add_argument("--species", required=True, help="Species name")
    build_parser.add_argument("--reference-dir", default=None, help="Reference directory")
    build_parser.add_argument("--force", action="store_true", help="Force rebuild")

    args = parser.parse_args()

    if args.command == "harmonize":
        from stangene import run
        try:
            result = run(
                path=args.input, species=args.species,
                output_dir=args.output_dir, dataset_name=args.dataset_name,
                reference_dir=args.reference_dir,
            )
            print(f"Harmonization complete. {len(result.mapping_table)} features processed.")
            print(f"Status counts: {result.stats}")
            print(f"Reports written to: {args.output_dir}")
        except Exception as e:
            logger.error("Harmonization failed: %s", e)
            sys.exit(1)
    elif args.command == "build-refs":
        from stangene.references import build_reference
        try:
            build_reference(species=args.species, reference_dir=args.reference_dir, force=args.force)
            print(f"Reference build complete for {args.species}.")
        except Exception as e:
            logger.error("Reference build failed: %s", e)
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
