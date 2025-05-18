parameters ={"input_seq":'/Users/devyaniravi/Downloads/result (12).fasta', 
             "output_folder":'/Users/devyaniravi/Desktop/toeholdoutputs',
             "reporter":'',
             "mol_type":'RNA',
             "min_unpaired":7,
             "suitable_aa":["G", "A", "S", "P", "V", "T", "C"]
            }   
from Bio import SeqIO
from Bio.Seq import Seq
import os
import numpy as np
import pandas as pd
import csv
import ViennaRNA
import matplotlib.pyplot as plt
import forgi
import forgi.visual.mplotlib as fvm
import forgi.graph.bulge_graph as cgb
class WrongFileExtension(Exception):
    """Exception raised for errors in the input file.

    Attributes:
        extension -- input file extension
        message -- explanation of the error
    """

    def __init__(self, extension, message="File must be in fasta format"):
        self.extension = extension
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.extension} -> {self.message}'
    
def DNAtoRNA(dnaseq):
    ''' Translate a DNA sequence to an RNA sequence by replacing T with U. 
    
    dnaseq (str): DNA sequence

    return : RNA sequences (str)    
    '''
    rnaseq=dnaseq.upper().replace('T', 'U')
    return rnaseq


def RNAtoDNA(rnaseq):
    ''' Translate an RNA sequence to a DNA sequence by replacing U with T. 

    rnaseq (str): RNA sequence

    return : DNA sequences (str)
    '''
    dnaseq=rnaseq.upper().replace('U', 'T')
    return dnaseq


def translate(rnaseq):
    ''' Translate an RNA sequence to an amino-acid sequence if there is a start codon. 
       
    rnaseq (str): RNA sequence 
    
    return :  amino-acid sequence (str) 
    '''
    
    # RNA codon table
    table = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
        "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
        "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
       "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
       "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
        "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
           "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
           "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
        }
    
    # Loop through the RNA sequence to find an AUG codon and store its position
    protein =""
    startpos="NaN"
    for i in range(0, len(rnaseq), 1):
        codon = rnaseq[i:i + 3]
        if codon=="AUG":
            startpos=i
            
        if startpos != "NaN":
            break
    
    # If a start codon is found, translate the RNA sequence to a protein using the AUG codon position and the RNA codon table 
    if startpos != "NaN":
        for i in range(startpos, len(rnaseq), 3):
            codon = rnaseq[i:i + 3]
            protein+= table[codon]
        
    return protein    
    
    

def get_switch_recognition_seq(trigger):
    ''' This function receives a target trigger miRNA sequence and obtains the recognition sequence for it.
    '''
    trigger_seq = Seq(trigger)
   
    return(trigger_seq.back_transcribe().reverse_complement().transcribe())


def write_toehold(seq_folder,toehold_seq,toehold_type,index):
    ''' Create a .txt file containing the toehold sequence.
    
    
    seq_folder (str): path where the sequence is contained
    
    toehold_seq (str): toehold sequence
    
    toehold_type (str): molecule type (DNA or RNA)
    
    index (int): toehold id 
   
    '''
        
    output_file=os.path.join(seq_folder, 'toehold_{!s}_{!s}.txt'.format(index,toehold_type))
    handle = open(output_file, 'w')
    handle.write('1\n')
    handle.write(str(toehold_seq) + '\n')
    handle.write('1')
    handle.close()


def get_full_switch(recognition_sequence_unpaired, recognition_sequence_paired, reporter, length_paired):
    ''' This sequence receives the recognition sequence of both the unpaired and the paired region of the toehold switch and the reporter gene to design the 
    full toehold switch. 
    
    recognition_sequence_unpaired : str

    recognition_sequence_paired : str
    
    reporter : str
    
    length_unpaired : int 
    
    
    return : str
    
    '''
    
    # Define the rest of the parts of the sequences
    part1 = 'GG'

    # The recognition sequence in the miRNA is divided in a unpaired and a paired part to the toehold structure
    part2_unpaired = recognition_sequence_unpaired
    part2_paired = recognition_sequence_paired
    
    # Get the reverse complement of the trigger to close the loop.
    rev_comp = str(part2_paired.reverse_complement())
    
    # Parts 3 and 5 should be complementary
    # Adjust their length so that the bottom stem has 11 total paired bases with 3 pairs at the bottom being 1 G/C and 2 A/U
    
    # Get the number of needed bases in the bottom stem
    needed = 11 - length_paired
    
    # Sequence before ribosome binding site in the top stem
    part3 = 'CUUUA'

    # Sequence of the ribosome binding site
    part4 = 'GAACAGAGGAGA'

    # Sequence after ribosome binding site in the top stem
    part5 = 'UAAAG'

    # Recognition site with added 'U' to make the bottom stem 11 with 'GGA' representing the sequence for wobble base pairing with the start codon
    part2_paired = part2_paired + needed * 'U' + 'GGA'

    # Start codon with added 'A' to make the bottom stemm 11 with the revese complementary sequence of the recognition site of the miRNA
    part6 = 'AUG' + needed * 'A' + rev_comp 
    

    # Linker sequence
    part7 = 'AAACCTGGCGGCAGCGCAAAAG'
    
    
    # Test for stop codons and start codon
    stop_codons = ['UGA', 'UAA', 'UAG','AUG'] 
    test_region = part6[3:] + part7
    
    # Split the test region into codons to check one at a time
    test_region_codons = [test_region[i:i+3] for i in range(0, len(test_region), 3)]
                          
    # Loop through the list and check if there are any stop codons in this reading frame
    for codon in test_region_codons:
        # assert not codon in stop_codons, 'The generated switch contains a stop codon' 
        if codon in stop_codons:
            return('Stop')
                          
    part8 = reporter
    
    full_toehold = part1 + part2_unpaired + part2_paired + part3 + part4 + part5 + part6 + part7 + part8
    
    return(full_toehold)

def testprotein(protein,nb,suitable_aa):
    ''' This function receives the protein sequence, the number of aa to be tested and the list with the aa wanted.
    
    protein (str): amino acid sequence in capital letter
    
    nb (int): number of amino acid to test at the beginning of the protein sequence 
    
    suitable_aa (list): amino acid in capital letter that are wanted at the beginning of the protein
    
    
    return : boolean (True : all amino acids tested were in the list. False : one amino acid tested was not in the list)
    
    '''
    
    suitable=True
    if protein[0]!="M":
    
        suitable=False
    elif nb > 0:
        for i in range(1,nb+1):
            if suitable==True:
                if protein[i] not in suitable_aa:
                    #print(f'Amino Acid {i+1} = {protein[i]} - AA not favorable \nSequence discarded \n ')
                    suitable=False
    return suitable
        

def SwitchMiDesigner_adjusted(parameters):
    ''' Main function to create the toehold switch. Necessite functions from SwitchDesigner_helper_functions.py.
    If not made yet, it creates a folder named after the miRNA. This folder will contain:
                .txt file containing the input variables, 
                .csv file for toehold switch candidates, 
                figures with the predicted secondary structures
                .txt files containing sequence (RNA or DNA) for each toehold switch.
    
    parameters (dict) ; must contain the following keys :
            "input_seq", "output_folder", "min_unpaired", "reporter", "mol_type", "suitable_aa"
    
    
    input_seq (str): path to the miRNA sequence (fasta file)
    
    output_folder (str): path to the output folder
    
    min_unpaired (int): Minimum number of unpaired residues in the secondary structure of the target miRNA for a candidate trigger to be considered
    
    reporter (str): DNA or RNA sequence - Reporter gene or tag to be added to the end of the toehold
    
    mol_type (str): 'DNA' or 'RNA'
    
    suitable_aa (list): Capital letters of amino acids that are wanted for the beginning of the protein sequence
    
    return : dataframe of selected toehold switch with all information collected.
    
    '''
    ### 1. Prepare input parameters
    # Retrieve parameters
    input_seq=parameters["input_seq"]
    output_folder=parameters["output_folder"]
    min_unpaired=parameters["min_unpaired"]
    reporter = parameters["reporter"]
    mol_type = parameters["mol_type"]
    suitable_aa = parameters["suitable_aa"]

    # Obtain the format of the input file containing the miRNA sequence
    extension=input_seq.split('.')[-1]

    # Raise exeption when input file is not a FASTA file
    if extension != "fasta":
        raise WrongFileExtension(extension)
        
    # Retrieve sequence name from the input_variable.py file
    name_sequence=input_seq.split('/')[-1][0:-6]

    # Create a folder "output" containing folders with each result. Create a second folder for each sequence
    seq_folder=output_folder+"/"+name_sequence
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(seq_folder):
        os.makedirs(seq_folder)
   
    # Make an input parameter file in the new directory
    f= open(seq_folder+"/input_variable.txt","w+")
    variable_str="Sequence : "+str(input_seq)
    variable_str=variable_str+"\nOutput_folder : "+str(output_folder)
    variable_str=variable_str+"\nMinimum unpaired residues : "+str(min_unpaired)
    variable_str=variable_str+"\nReporter : "+str(reporter)
    variable_str=variable_str+"\nmol_type : "+str(mol_type)
    variable_str=variable_str+"\nSuitable amino acids : "+str(suitable_aa)
    f.write(variable_str)
    f.close() 

    
    ### 2. Prepare the miRNA sequence
    # Parse the input sequence
    input_seq_record = next(SeqIO.parse(input_seq, 'fasta'))

    # If input sequence is DNA, parse and transcribe
    if mol_type == 'DNA':   
        seq_miRNA = DNAtoRNA(str(input_seq_record.seq))
        
    # If it is RNA, make sure it is in upper letters.
    elif mol_type == 'RNA':
        seq_miRNA = str(input_seq_record.seq).upper()

    # Calculate the secondary structure and MFE of miRNA with ViennaRNA
    structure_miRNA, mfe_miRNA = ViennaRNA.fold(seq_miRNA)
    print("ViennaRNA created and MFE proxy structure calculated")
    
    # Store values
    miRNA_dict={}
    miRNA_dict['sequence']=seq_miRNA
    miRNA_dict['structure']=structure_miRNA
    miRNA_dict['energy']=mfe_miRNA
        
    # Create list for toehold results
    results = []

    # Set index
    index=0
    
    print("miRNA Scanner to find trigger sequence")

    # Get the recognition sequence of the miRNA to implement in the toehold switch
    recognition_seq_miRNA = get_switch_recognition_seq(seq_miRNA)

    # Minimum length of miRNA that should bind in the unpaired region (before the toehold structure)
    length_unpaired = 4
    
    # Minumum length of miRNA that should bind in the paired region (in the toehold structure)
    length_paired = 3

     
    ### 3. Find candidates : Loop through the miRNA sequence to find a suitable hairpin start position
    for hairpin_start_pos in range(length_unpaired, len(miRNA_dict['sequence'])+1 - length_paired):
        index=index+1
        
        # Store the sequence of possible hairpin's base
        hairpin_base = recognition_seq_miRNA[hairpin_start_pos:hairpin_start_pos + 3]   
        
        # Discard cases for which the base of the hairpin does not have two weak bases and a strong base
        if hairpin_base.count('A') + hairpin_base.count('U') != 2:
            continue


        ### 4. Generate toehold 
        # Split recognition sequence in unpaired and paired region in the toehold switch
        unpaired_region =  recognition_seq_miRNA[0:hairpin_start_pos]
        paired_region =  recognition_seq_miRNA[hairpin_start_pos:len(seq_miRNA)]

        # Check if the number of nucleotides in miRNA paired region exceeds 11
        if len(paired_region) <= 11:
            # Paired region stays the same if smaller than 11
            paired_region =  paired_region

            # State that there will be no nculeotides that are unable to bind to the toehold switch
            ntmiRnotbound = 'none'

            # Give the miRNA sequence and structure that will be able to bind to the toehold
            seq_bound_miRNA = seq_miRNA
            structure_bound_miRNA = structure_miRNA
            
        elif len(paired_region) > 11:
            # Paired region in miRNA will be max. 11 nucleotides
            end_position = len(unpaired_region) + 11

            # State how many nucleotides are unable to bind to the toehold switch
            ntmiRnotbound = len(seq_miRNA) - end_position

            # Obtain sequence of paired region of max. 11 nucleotides
            paired_region = recognition_seq_miRNA[hairpin_start_pos:end_position]
            
            # Give the miRNA sequence and structure that will be able to bind to the toehold
            seq_bound_miRNA = seq_miRNA[ntniRnotbound:len(seq_miRNA)]
            structure_bound_miRNA = structure_miRNA[ntniRnotbound:len(seq_miRNA)]

        
        # Count the number of unpaired bases in the miRNA able to anneal to the toehold switch
        notannealed = structure_bound_miRNA.count('.')
        
        # Discard cases for which the number of unpaired bases in the subsequence doesn't match the minimum number of unpaired bases.
        if notannealed < min_unpaired:
            #print('not enought unpaired bases in the subsequence')
            continue          

        # Obtain length of paired region
        length_paired = len(paired_region)

        # Generate toehold switch
        seq_toehold = get_full_switch(unpaired_region, paired_region, reporter, length_paired)
            
        # Store values if toehold is generated
        toehold_dict={}
        if seq_toehold != 'Stop':
            seq_toehold=str(seq_toehold)
            seq_toehold=DNAtoRNA(seq_toehold)
            seq_toehold_DNA=RNAtoDNA(seq_toehold)
        
            # Calculate the secondary structure and MFE of toehold switch with ViennaRNA
            structure_toehold, mfe_toehold = ViennaRNA.fold(seq_toehold)

            # Store values of toehold switch in dictionary
            toehold_dict["sequence_toehold"] = seq_toehold
            toehold_dict['DNA_sequence_toehold'] = seq_toehold_DNA
            toehold_dict['structure_toehold'] = structure_toehold
            toehold_dict['mfe_toehold'] = mfe_toehold

        # Check and save properties of toehold switch
        if toehold_dict:
            # Store values of miRNA sequence and index
            toehold_dict["input_miRNA"]=seq_miRNA
            toehold_dict["index"] = index

            ### 5. Check how well the toeholds bind to the miRNA
            # Hybridize miRNA to the toehold switch
            hybrid_sequence = seq_toehold + '&' + seq_miRNA

            # Calculate the secondary structure and MFE of toehold switch with ViennaRNA
            structure_hybrid, mfe_hybrid = ViennaRNA.cofold(hybrid_sequence)
            
            # Discard toehold if binding MFE is bigger than the sum of each when not bound
            if mfe_hybrid > (mfe_toehold + mfe_miRNA):
                continue
            
            # Calculate the GC content of the toehold switch
            GC_content = round(float(seq_toehold.count('G') + seq_toehold.count('C'))*100/len(seq_toehold), 2)

            ### 6. Verify the protein produced by the toehold switch
            # Translate the protein
            protein=translate(toehold_dict["sequence_toehold"])
            
            # Discard cases for which the protein generated by the reporter 1. do not begin by M and 2. the 2nd and 3rd amino acids are not amino acids with low molecular weight.
            suitable=testprotein(protein,3,suitable_aa)
            if suitable==False:
                aminoacids = "amino_acids_not_from_input"
            else:
                aminoacids = "amino_acids_from_input"


            ### 7. Generate figures of the predicted secondary structure
            # Parse the sequence and structure into a forgi BulgeGraph object
            bg = cgb.BulgeGraph.from_dotbracket(structure_toehold, seq_toehold)

            # Create the plot of secondary structure
            fig, ax = plt.subplots(figsize=(10, 10))
            fvm.plot_rna(bg, text_kwargs={"fontsize": 12}, ax=ax)
            plt.tight_layout()

            # Save the plot as PNG
            toehold_png_file = os.path.join(seq_folder, f"toehold_{index}_ViennaRNA.png")
            plt.savefig(toehold_png_file)
            plt.close()

            # Save the RNA and DNA sequence in a text file  
            write_toehold(seq_folder,toehold_dict["sequence_toehold"],'RNA',index)
            write_toehold(seq_folder,toehold_dict['DNA_sequence_toehold'],'DNA',index)

            # Expand the results list
            results_list = [
                index,
                seq_miRNA,           
                structure_miRNA, 
                mfe_miRNA,
                seq_bound_miRNA,
                ntmiRnotbound,
                structure_bound_miRNA,
                notannealed,
                seq_toehold,
                structure_toehold, 
                mfe_toehold,
                structure_hybrid, 
                mfe_hybrid,
                GC_content,
                protein,
                aminoacids]
            results.append(results_list)

    
    ### 7. Create csv files for results
    if len(results)>0:
        results_df = pd.DataFrame(results, columns=['Index', 'Sequence_miRNA', 'Structure_miRNA',  'MFE_miRNA(kcal/mol)', 'miRNA_Sequence_Binding_to_Toehold', 'Number_Nucleotides_Not_Bound_to_Toehold', 'Structure_miRNA_Bound_to_Toehold', 'Number_Unpaired_Bases_miRNA_Subsequence', 'Sequence_Toehold', 'Structure_Toehold', 'MFE_Toehold_Switch(kcal/mol)', 'Structure_Hybrid', 'MFE_Hybrid(kcal/mol)', 'GC_Content_Toehold', 'Protein', 'Amino_Acids_from_Input'])
        sorted_results = results_df.sort_values(['MFE_Toehold_Switch(kcal/mol)'], ascending = False)
        sorted_results.to_csv(path_or_buf=os.path.join(seq_folder, 'selected_toeholds_results.csv'), sep = '\t', index = False)
        print("Toehold canditates stored of ViennaRNA")
    
    else: 
        print("No toehold canditates")

SwitchMiDesigner_adjusted(parameters)
