### Implementing and Evaluating Byte Pair Encoding (BPE) for DNA Tokenization

############ Project 7 ############ 

filename="/home/ioanna/metaptixiako/kanterakis/chr20.fa.gz"

####Preprocessing steps
import gzip

def preprocess_file(filename):
    sequence_parts=[]
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                continue
            clean_line=line.strip()
            clean_upper_line=clean_line.upper()
            sequence_parts.append(clean_upper_line)
        full_sequence="".join(sequence_parts)
    
    return full_sequence


####Part 1
def count_pairs(tokens): #from a list of tokens it returns a dict with the count of each pair
    pairs={}
    for i in range(len(tokens)-1):
        pair=(tokens[i],tokens[i+1])
        if pair in pairs:
            pairs[pair]+=1
        else:
            pairs[pair]=1
    
    return pairs


def merge_pairs(tokens, pair_to_merge, new_token): #it takes a list of tokens, the pair of tokens that we want to merge in one and the new token that will replace them. it returns a list of the new tokens
    new_tokens=[]
    i=0
    while i<len(tokens):
        if i==len(tokens)-1: #if we are in the last possition then we just add this token
            new_tokens.append(tokens[i])
            i+=1
            continue
        if (tokens[i],tokens[i+1])==pair_to_merge: #we check if two tokens that are next to each other are the pair we want to merge
            new_tokens.append(new_token) #if it is then we put the new token on the list 
            i+=2 #we skip these two tokens we just checked
        else:
            new_tokens.append(tokens[i]) #if it is not we put the current token in the list
            i+=1 #we skip this one token
    
    return new_tokens


def train_bpe(sequence, K):
    tokens=list(sequence) #the list of tokens
    final_tokens=set(tokens) #the set of the final tokens (it already contains the individual characters/bases as tokens)  
    
    initial_tokens=len(final_tokens)
    if K<=initial_tokens: #if K<=from initial tokens, then we return the initial tokens
        return final_tokens 

    rest_tokens=K-initial_tokens #how many more tokens we need to find
    
    for i in range(rest_tokens):
        counts=count_pairs(tokens) #we count the pairs
        if not counts: #if we have only one token, then we stop
            break
        pair_to_merge=max(counts, key=counts.get) #the pair that we will merge is the one with the biggest count
        new_token=pair_to_merge[0]+pair_to_merge[1] #we make the new token
        final_tokens.add(new_token) #we add the new token to the list of the final tokens
        tokens=merge_pairs(tokens, pair_to_merge, new_token) #we update the list of tokens with the merged pair
    
    return final_tokens



####Part 2
def tokenize_dna(sequence, bpe_vocab): #it takes a sequence and a set of tokens. it returns a list of tokens that represent the sequence
    tokenized_sequence=[]
    sorted_vocab=sorted(list(bpe_vocab), key=len, reverse=True) #it sorts the tokens from biggest to smallest (we want to have the biggest pair possible)

    i=0
    while i<len(sequence):
        found_match=False
        for token in sorted_vocab: #we try each token from the sorted list of tokens
            if sequence[i:].startswith(token): #if the sequence in the position i starts with the token
                tokenized_sequence.append(token) #we add the token in the list
                i+=len(token) #we move after the bases that make up the token
                found_match=True
                break #we stop the for loop and continue the while loop for the next position of the sequence
        if not found_match: #if we do not find a match we add the base in this position and continue on the next 
            tokenized_sequence.append(sequence[i])
            i+=1
    return tokenized_sequence



####Part 3
import random

def evaluate_tokenization(sequence, K, M, runs): #it returns the average number of tokens per sequence
    max_start_index=len(sequence)-M #the max idex position to get a sequence of length M
    if max_start_index<0:
        print("sub-sequence is bigger than the original sequence")

    total_tokens=0
    for i in range(runs):
        train_start=random.randint(0,max_start_index)
        train_sequence=sequence[train_start:train_start+M] #this is the sequence to train the bpe

        bpe_vocab=train_bpe(train_sequence, K) #this is the vocab we get after we train the bpe

        test_start=random.randint(0,max_start_index)
        while test_start==train_start:
            test_start=random.randint(0,max_start_index) #we do not want to get the same test sequence that we trained the bpe with
        test_sequence=sequence[test_start:test_start+M] #this is the sequence to test the bpe

        tokenized_sequence=tokenize_dna(test_sequence, bpe_vocab) #tokenization of the sequence
        total_tokens+=len(tokenized_sequence) #we add the count of the tokens 
    
    average_tokens=total_tokens/runs #this is the average number of tokens per sequence
    return average_tokens



####Part 4
parameter_k=[10,50,100,200,500,1000]
parameter_m=[1000,5000,10000,50000]
sequence=preprocess_file(filename)
runs=100

results=[]

for k in parameter_k:
    for m in parameter_m:
        print(f'Evaluating for K={k} and M={m}')
        avg=evaluate_tokenization(sequence,k,m,runs)
        results.append({'K':k, 'M':m, 'Average Tokens':avg})
        print(f'Average tokens={avg}')

print('K,   M,  Average Tokens')
for res in results:
    print(res['K'],   res['M'], res['Average Tokens'])