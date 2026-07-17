import csv
import subprocess
import click

# Path to your script that takes multiple inputs
script_path = "/broad/macosko/mkim/Scripts/230720_Synapseseq_Slideseq_Dialout_Filtering.py"

@click.command()
@click.argument("input-csv", type=click.Path(exists=True))

def main(
    input_csv,
):
    with open(input_csv, "r") as file:
        # Create a CSV reader object
        reader = csv.DictReader(file)

        # Iterate over each row in the CSV file

        i=1
        for row in reader:

            print(f"Processing sample #{i}")

            # Build the command to run your script with the extracted value
            print(f"Setting up job")
            command = ["python", script_path, \
                row['input_path'], \
                "--input-sample", row['input_sample'], \
                "--output-folder", row['output_folder'], \
                "--polygon-folder", row['polygon_folder'], \
                "--polygon-list", row['polygon_list']
            ]

            print(f"Submitting job")
            # Run the command
            subprocess.run(command)
            print(f"Job completed")
            i=i+1

if __name__ == "__main__":
    main()