import random
import numpy as np

# Define constants
NUM_CHROMOSOMES = 23
CHROMOSOME_LENGTH = 10000
NUM_GENES = 20
MEAN_GENE_LENGTH = 100
STD_GENE_LENGTH = 500
MUTATION_RATE = 0.1

# Define probability distributions
GC_CONTENT_MEAN = 0.7
GC_CONTENT_STD = 0.1
GENE_ORIENTATION_PROB = 0.5
GENE_OVERLAP_PROB = 0.1
PSEUDOGENE_PROB = 0.05
REPETITIVE_ELEMENT_PROB = 0.5

# Define regulatory elements
REGULATORY_ELEMENT_TYPES = ['promoter', 'enhancer', 'silencer', 'insulator', 'repressor', 'activator', 'cofactor']
REGULATORY_ELEMENT_LENGTH = 500
specific_mutations = [
    {'position': 50, 'base': 'G'},  # replace A with G at position 500
    {'position': 10, 'base': 'T'},  # replace C with T at position 1000
    {'position': 20, 'base': 'A'}  # replace G with A at position 2000
]
def generate_dna_sequence(length, gc_content_mean, gc_content_std):
    """Generate a DNA sequence of a specified length with a specified GC content."""
    gc_content = np.random.normal(gc_content_mean, gc_content_std)
    gc_content = max(0, min(1, gc_content))  # Ensure GC content is between 0 and 1
    dna_sequence = ''
    for _ in range(length):
        if random.random() < gc_content:
            dna_sequence += random.choice('GC')
        else:
            dna_sequence += random.choice('AT')
    return dna_sequence

def generate_gene(chrom_sequence, start, end):
    """Extract a gene sequence from a chromosome sequence."""
    gene_sequence = chrom_sequence[start:end]
    return gene_sequence, start, end

def generate_regulatory_elements(chrom_sequence, gene_start, gene_end):
    """Generate regulatory elements (e.g., promoters, enhancers, silencers) near a gene."""
    regulatory_elements = []
    for _ in range(random.randint(1, 3)):
        element_type = random.choice(REGULATORY_ELEMENT_TYPES)
        element_length = random.randint(10, 50)
        element_start = random.randint(gene_start - 100, gene_start - 1)
        element_end = element_start + element_length
        if element_end > gene_start:
                   continue
        element_sequence = chrom_sequence[element_start:element_end]
        regulatory_elements.append((element_type, element_start, element_end))
    return regulatory_elements

def generate_repetitive_elements(chrom_sequence, gene_start, gene_end):
    """Generate repetitive elements (e.g., Alu, LINE) near a gene."""
    repetitive_elements = []
    for _ in range(random.randint(1, 3)):
        element_length = random.randint(10, 50)
        element_sequence = ''.join(random.choice('ATGC') for _ in range(element_length))
        element_start = random.randint(gene_start - 100, gene_end)
        repetitive_elements.append((element_sequence, element_start))
    return repetitive_elements

# Generate chromosome sequences
chrom_sequences = [generate_dna_sequence(CHROMOSOME_LENGTH, GC_CONTENT_MEAN, GC_CONTENT_STD) for _ in range(NUM_CHROMOSOMES)]

# Open a file to write the genome sequence
with open('synthetic_genome.fasta', 'w') as f:
    for i, chrom_sequence in enumerate(chrom_sequences):
        print(f"Writing chromosome {i+1} to file...")
        f.write(f'>Chromosome {i + 1}\n')
        f.write(f'{chrom_sequence}\n')

        gene_starts = list(range(0, len(chrom_sequence), CHROMOSOME_LENGTH // NUM_GENES))
        random.shuffle(gene_starts)

        for gene_start in gene_starts:
            gene_end = min(gene_start + CHROMOSOME_LENGTH // NUM_GENES, len(chrom_sequence))
            gene_sequence, start, end = generate_gene(chrom_sequence, gene_start, gene_end)

            print(f"Writing gene {gene_start} to file...")
            f.write(f'>Forward gene {gene_start}\n')
            f.write(f'{gene_sequence}\n')

            # Introduce specific mutations
            for mutation in specific_mutations:
                if mutation['position'] >= start and mutation['position'] < end:
                    gene_sequence = gene_sequence[:mutation['position'] - start] + mutation['base'] + gene_sequence[mutation['position'] - start + 1:]
                    f.write(f'>{mutation["base"]} mutation at position {mutation["position"]}\n')
                    f.write(f'{gene_sequence}\n')

            regulatory_elements = generate_regulatory_elements(chrom_sequence, gene_start, gene_end)
            for element_type, element_start, element_end in regulatory_elements:
                print(f"Writing regulatory element {element_type} {element_start} to file...")
                f.write(f'>{element_type} {element_start}\n')
                f.write(f'{chrom_sequence[element_start:element_end]}\n')

            repetitive_elements = generate_repetitive_elements(chrom_sequence, gene_start, gene_end)
            for element_sequence, element_start in repetitive_elements:
                print(f"Writing repetitive element {element_start} to file...")
                f.write(f'>Repetitive element {element_start}\n')
                f.write(f'{element_sequence}\n')

                # Add pseudogenes
                if random.random() < PSEUDOGENE_PROB:
                   pseudogene_sequence = ''.join(random.choice('ACGT') for _ in range(MEAN_GENE_LENGTH))
                   f.write(f">Pseudogene {element_start+1}\n")
                   f.write(pseudogene_sequence + "\n")

                # Add repetitive elements
                if random.random() < REPETITIVE_ELEMENT_PROB:
                   repetitive_element_sequence = ''.join(random.choice('ACGT') for _ in range(MEAN_GENE_LENGTH))
                   f.write(f">Repetitive element {element_start+1}\n")
                   f.write(repetitive_element_sequence + "\n")

            # Add regulatory elements to the gene
            for k in range(random.randint(1, 5)):  # Add 1-5 regulatory elements per gene
                regulatory_element_type = random.choice(REGULATORY_ELEMENT_TYPES)
                regulatory_element_sequence = ''.join(random.choice('ACGT') for _ in range(REGULATORY_ELEMENT_LENGTH))
                f.write(f">{regulatory_element_type} {k+1}\n")
                f.write(regulatory_element_sequence + "\n")

            # Add additional regulatory elements for better reinforcement
            for k in range(random.randint(1, 3)):  # Add 1-3 additional regulatory elements per gene
                additional_regulatory_element_type = random.choice(['enhancer', 'silencer', 'cofactor'])
                additional_regulatory_element_sequence = ''.join(random.choice('ACGT') for _ in range(REGULATORY_ELEMENT_LENGTH))
                f.write(f">{additional_regulatory_element_type} {k+1}\n")
                f.write(additional_regulatory_element_sequence + "\n")

            # Write the gene sequence
            f.write(gene_sequence + "\n")

print("Synthetic genome sequence written to synthetic_genome.fasta")
           