import nupack

# Define the mRNA sequence
rna_sequence = "AUGACCGAACAGCAGUGGAACUUUGCGGGCAUUGAAGCGGCGGCGAGCGCGAUUCAGGGCAACGUGACCAGCAUUCAUAGCCUGCUGGAUGAAGGCAAACAGAGCCUGACCAAACUGGCGGCGGCGUGGGGCGGCAGCGGCAGCGAAGCGUAUCAGGGCGUGCAGCAGAAAUGGGAUGCGACCGCGACCGAACUGAACAACGCGCUGCAGAACCUGGCGCGCACCAUUAGCGAAGCGGGCCAGGCGAUGGCGAGCACCGAAGGCAACGUGACCGGCAUGUUUGCG"

# Initialize the NUPACK model for RNA
model = nupack.Model(material='rna')

# Predict the minimum free energy (MFE) secondary structure with the model
structures = nupack.mfe(rna_sequence, model)

# Output the predicted structure and its MFE for each possible structure
for structure in structures:
    print("Predicted secondary structure:", structure.structure)
    print("Minimum Free Energy:", structure.energy)
