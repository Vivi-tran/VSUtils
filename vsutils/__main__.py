import argparse
import pandas as pd
from vsutils.ADMET.mc_filter import MCFilter
import os


def validate_input(input_arg):
    """Check if the input is a file path or a SMILES string."""
    if os.path.isfile(input_arg):
        return "csv"
    else:
        return "smiles"


def process_smiles(smiles, filter, output_csv=None):
    """Process a single SMILES string using MCFilter."""
    result = filter.calculate_des(smiles)
    result_df = pd.DataFrame([result])  # Convert Series to DataFrame for uniformity
    if output_csv:
        result_df.to_csv(output_csv, index=False)
        print(f"Results saved to {output_csv}")
    else:
        print("Results for the provided SMILES string:")
        print(result_df)


def process_csv(input_csv, smiles_column, filter, output_csv=None):
    """Process a CSV file containing SMILES strings using MCFilter."""
    df = pd.read_csv(input_csv)
    result_df = filter.process_dataframe(df, smiles_column)
    if output_csv:
        result_df.to_csv(output_csv, index=False)
        print(f"Processing complete. Results saved to {output_csv}")
    else:
        print("Results for the processed CSV data:")
        print(result_df)


def main():
    parser = argparse.ArgumentParser(
        description="Process SMILES strings or CSV files containing SMILES using MCFilter."
    )
    parser.add_argument("input", type=str, help="Input SMILES string or CSV file path.")
    parser.add_argument(
        "-c",
        "--smiles_column",
        type=str,
        default="smiles",
        help="Column name in CSV that contains the SMILES strings.",
        dest="column",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Optional output CSV file path to save the results.",
        dest="output_csv",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel jobs to run.",
        dest="n_jobs",
    )
    parser.add_argument(
        "-v", "--verbose", type=int, default=0, help="Verbosity level.", dest="verbose"
    )

    args = parser.parse_args()

    # Create MCFilter instance
    filter = MCFilter(n_jobs=args.n_jobs, verbose=args.verbose)

    input_type = validate_input(args.input)

    if input_type == "csv":
        process_csv(args.input, args.column, filter, args.output_csv)
    elif input_type == "smiles":
        process_smiles(args.input, filter, args.output_csv)


if __name__ == "__main__":
    main()
