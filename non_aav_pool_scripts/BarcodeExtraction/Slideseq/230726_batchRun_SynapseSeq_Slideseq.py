import csv
import subprocess
import click

# Path to your script that takes multiple inputs
script_path = "/broad/macosko/mkim/Scripts/230726_SynapseSeq_Slideseq_Dialout_NoDegen.py"

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

            # Make output directory
            output_dir = row['output_dir']
            subprocess.run(["mkdir", output_dir])
            print(f"Directory created at: {output_dir}")

            # Build the command to run your script with the extracted value
            print(f"Setting up job")
            command = ["python", script_path, \
                row['read1_path'], row['read2_path'], \
                "--puck-dir", row['puck_path'], \
                "--output-dir", output_dir, \
                "--tag-mismatch", row['tag_mm'], \
                "--const-mismatch", row['constant_mm'], \
                "--debug"
            ]

            print(f"Submitting job")
            # Run the command
            subprocess.run(command)
            print(f"Job completed")
            i=i+1

if __name__ == "__main__":
    main()