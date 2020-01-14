# bioinf-projekt
Project as part of the [Bioinformatics](https://www.fer.unizg.hr/predmet/bio) class at the [Faculty of Electrical Engineering and Computing](https://www.fer.unizg.hr/) at the [University of Zagreb] (http://www.unizg.hr/).

Authors: Andrej Mijić, Ilija Domislović, Zorica Žitko.

# Installation instructions
Run make.

# Running the program
Run mutation_checker with all 7 input arguments.

Ulazni argumenti:<br/>
kmer_size - Length of substring which will be indexed.<br/>
minimizer_size - Length of each minimizer.<br/>
number_of_threads - The number of threads the program will use.<br/>
mutation_voting_threshold - Threshold that defines how many aligned sequences need to vote for a mutation for that mutation to be entered into the final solution.<br/>
reference_genome_path - Path to reference genome file.<br/>
sequencing_results_path - Path to sequencing results file.<br/>
output_path - Path to output .csv file.<br/>
