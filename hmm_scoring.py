import subprocess
import pandas as pd
from Bio import SeqIO

#this function builds a hmm based on a Multiple sequences alignement. Takes an MSA.fastab file as first argument and the name of the desired output as second argument
def build_hmm(input_file, output_file):
    hmmbuild_command = f"hmmbuild {output_file} {input_file}"
    subprocess.call(hmmbuild_command, shell=True) #The hmm is built with the hmmer program, and this command calls the program via the terminal.
    
#This function computes the score of a sequence file against an hmm model. It takes a fasta file containing as many sequences as we want as first argument and an hmm model as second argument and returns the score, evalue and dame of each sequence as lists.
  def score_sequence(sequence_file, hmm_file):
    evalues = {}
    sequences = {}
    
    result = subprocess.run(['hmmsearch', hmm_file, sequence_file],text = True, capture_output= True, check=True) #The hmmsearch command from hmmer is used to compute the scores

    output = result.stdout
    lines = output.split("  ------- ------ -----    ------- ------ -----   ---- --  --------         -----------")[1].split(" ------ inclusion threshold ------")[0].split("\n") #This line processes the output of the hmmesearch command to onlyy keep the results matrix
    lines = lines[1:len(lines)-1]
    evalues = []
    scores = []
    sequences = []
    
    for line in lines: #this loop extracts the features we want from the results matrix
        evalues.append(float(line.split()[0]))
        sequences.append(line.split()[-1])
        scores.append(float(line.split()[1]))
    return evalues, sequences, scores 
  
def build_scoring_matrix(sequences_files, MSA_files, null_value = 0.5): # This function buils a scoring matrix that will be used to train our classification model. 
  #It takes as input a list of sequences.fasta files. Each file in the list corresponds to one of the four groups that the classification will try to differenciate. The second 
  #input is a list of MSA files, each MSA is based on one of the sequences files from the first output. The two lists have to be in the same order. The last argument is the value that will replace null values in our final matrix to avoid bias when training the model.
    
  scoring_matrix = pd.DataFrame({})
    HMM_files = []
    for index, file in enumerate(MSA_files): # The first loop calls the build_Hmm function to build an hmm by MSA file
        HMM_files.append(MSA_files[index].split(".")[0] + ".hmm")
        build_hmm(MSA_files[index], HMM_files[index])
    
    for hmm_index, hmm_file in enumerate(HMM_files): # Thes two nested loops compute the scores and e value for each sequences file against each hmm, and also stores the sequences names, their names and the group they belong to (ground truth to train the model on)
        all_evalues = []
        all_scores = []
        row_names = []
        labels = []
        lengths = []
        for seq_index, seq_file in enumerate(sequences_files):
            evalues, seq_names, scores = score_sequence(seq_file, hmm_file)
            all_evalues.extend(evalues)
            row_names.extend(seq_names)
            all_scores.extend(scores)   
            ground_truth = [seq_file.split(".")[0]] * len(evalues)
            labels.extend(ground_truth)
            if hmm_index == len(HMM_files) - 1:
                for single_seq in SeqIO.parse(seq_file, "fasta"):
                    if single_seq.id in row_names:
                        lengths.append(len(single_seq.seq))     
       
        scoring_matrix[hmm_file.split(".")[0] + str(hmm_index) + " Evalue"] = all_evalues
        scoring_matrix[hmm_file.split(".")[0] + str(hmm_index) + " score"] = all_scores
    scoring_matrix["lengths"] = lengths
    scoring_matrix["labels"] = labels
    scoring_matrix.index = row_names
    scoring_matrix.replace(0, null_value, inplace=True)
    
    return scoring_matrix

# Here we will train a model to classify SNARE sequences in four groups "Qa, Qb, Qc and R". 
matrix = build_scoring_matrix(["archeaplastida_SNARE.fasta","archeaplastida_SNARE.fasta","archeaplastida_SNARE.fasta","archeaplastida_SNARE.fasta"],['Output_AViri_QA.fasta','Output_AViri_QA.fasta','Output_AViri_QA.fasta','Output_AViri_QA.fasta'], 0.5)
matrix.to_csv('matrix_scoring.csv')
