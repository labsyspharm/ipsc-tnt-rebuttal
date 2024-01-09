import click
import pandas as pd
import os
import shutil

@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.argument('accession_column', type=str)
def merge_fastq(input, accession_column):
    # Load the TSV file into a pandas DataFrame
    df = pd.read_csv(input, sep='\t')

    # Group by sample_accession
    groups = df.groupby(accession_column)

    # Iterate over each group
    for sample_accession, group_df in groups:
        # Create new files for the merged data
        output_file_1 = f'{sample_accession}_1.fastq.gz'
        output_file_2 = f'{sample_accession}_2.fastq.gz'
        output_file = f'{sample_accession}_1.fastq.gz'

        click.echo(f'Merging files for sample accession {sample_accession}...')

        # Check if the output file already exists, and empty it if it does
        for file in [output_file_1, output_file_2, output_file]:
            if os.path.exists(file):
                with open(file, 'wb'):
                    pass
                click.echo(f'Emptied existing file {file}')

        # Iterate over each row in the group
        for index, row in group_df.iterrows():
            run_accession = row['run_accession']
            input_file_1 = f'{run_accession}_1.fastq.gz'
            input_file_2 = f'{run_accession}_2.fastq.gz'
            input_file = f'{run_accession}.fastq.gz'

            # Merge the input files into the output files if they exist
            for input_file, output_file in zip([input_file_1, input_file_2, input_file], [output_file_1, output_file_2, output_file]):
                try:
                    with open(input_file, 'rb') as f_in, open(output_file, 'ab') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    click.echo(f'Merged {input_file} into {output_file}')
                except FileNotFoundError:
                    # Handle the case where the file doesn't exist
                    click.echo(f'File {input_file} not found, skipping...')

if __name__ == '__main__':
    merge_fastq()
