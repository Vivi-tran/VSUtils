# VSUtils

## Introduction
`VSUtils` is a Python-based utility designed for facilitating virtual screening processes in the realm of drug discovery. This tool incorporates `MCFilter`, a powerful module that offers extensive molecular property calculations and filtering capabilities based on established pharmaceutical rules such as Ro5, Pfizer, and GSK, among others.

## Installation
To set up `VSUtils` for use, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/VSUtils.git
   cd VSUtils
   ```

2. Install required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

3. Install vsutils packages:
   ```bash
   pip install vsutils
   ```
## Usage

### In a Jupyter Notebook
For interactive use, such as in a Jupyter notebook, follow this example to utilize the `MCFilter`:

```python
from vsutils.ADMET.mc_filter import MCFilter

# Initialize the MCFilter with default settings
filter = MCFilter()

# Calculate descriptors for a single SMILES string
smiles = "CCO"
results = filter.calculate_des(smiles)
print(results)
```

### Command-Line Interface
`MCFilter` can also be executed from the command line to process either individual SMILES strings or batches from a CSV file:

- **Processing a Single SMILES String:**
  ```bash
  python main.py "CCO"
  ```

- **Processing SMILES from a CSV File:**
  ```bash
  python main.py path/to/input.csv --smiles_column column_name --output path/to/output.csv
  ```

#### Command-Line Options
- `input`: Input SMILES string or path to a CSV file containing SMILES strings.
- `-c, --smiles_column`: Specifies the column name in the CSV that contains the SMILES strings. Default is 'smiles'.
- `-o, --output`: Specifies the path for the output CSV file. If not provided, results will be printed to the console.
- `-j, --jobs`: Number of parallel jobs to run for processing. Default is 1.
- `-v, --verbose`: Verbosity level of output. Default is 0 (silent).

## Requirements
This utility requires Python 3.9 or later. All dependencies are listed in `requirements.txt`.

## License
`VSUtils` is open-source software licensed under the MIT License - see the [License](LICENSE).

