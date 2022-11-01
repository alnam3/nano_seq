"""
  Author: Kristoffer Sahlin
"""
import os,sys
import random
import itertools
import argparse
import errno
import math


def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

'''
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def simulate_reads( args, recombinants ):

    reads = {}
    if args.X == 4:
        error_lvls = [0.9, 0.95, 0.96, 0.98, 0.99, 0.995]
    elif  args.X == 12:
        error_lvls = [0.75, 0.85, 0.875, 0.91, 0.95, 0.98]        
    elif args.X == 7:
        error_lvls = [0.85, 0.875, 0.9, 0.92, 0.96, 0.98, 0.99, 0.995]
    else:
        print("Wrong --X (error rate specified), choose between 4, 7 or 12 (%), or contact Kristoffer for other level.")
        sys.exit()

    pool_of_rec = list(recombinants.items())
    for i in range(args.N):
        rec_acc, recombinant = random.choice(pool_of_rec)
        read = []
        qual = []

        was_del = False
        for l, n in enumerate(recombinant):
            p_correct_reading = random.choice(error_lvls)
            p_error = 1.0 - p_correct_reading

            r = random.uniform(0,1)
            if r > p_correct_reading:
                error = True
            else:
                error = False

            if error:
                r = random.uniform(0,1)
                if r < 0.6: #deletion
                    was_del = p_error
                    pass 
                elif 0.6 <= r < 0.9:
                    read.append(random.choice("ACGT"))
                    qual.append( round(-math.log(p_error,10)*10) )

                else:
                    read.append(n)
                    qual.append( round(-math.log(p_error,10)*10) )

                    r_ins = random.uniform(0,1)
                    while r_ins >= 0.7:
                        read.append(random.choice("ACGT"))
                        r_ins = random.uniform(0,1)
                        qual.append( round(-math.log(0.7,10)*10) )

            else:
                if was_del: # adding uncertainty from prevous deleted base
                    read.append(n)
                    qual.append( round(-math.log(was_del,10)*10) )
                else:
                    read.append(n)
                    qual.append( round(-math.log(p_error,10)*10) )
                was_del = False

        if not read:
            continue
        read_seq = "".join([n for n in read])
        qual_seq = "".join([chr(q + 33) for q in qual])
        reads[str(rec_acc) + "_readID_" + str(i)  ] = (read_seq, qual_seq)

    return reads

    # if is_fastq:
    #     for acc, (read_seq,qual_seq) in sorted(reads.items(), key = lambda x: len(x[1]), reverse = True):
    #         outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, read_seq, "+", qual_seq))
    # else:
    #     for acc, (read_seq,qual_seq) in sorted(reads.items(), key = lambda x: len(x[1]), reverse = True):
    #         outfile.write(">{0}\n{1}\n".format(acc, read_seq))
    
    # outfile.close()


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def simulate_rec(args, orig_rec):
    new_rec = orig_rec.copy()
    cnt = 0
    while len(new_rec) < args.C:
        rec_acc = random.choice(list(new_rec.keys()))
        # mutate
        rec_seq = new_rec[rec_acc]
        L = len(rec_seq)
        nr_muts = sum([1 if random.random() < args.S else 0 for i in range(L) ]) 
        if args.exp_type == 'subs':
            muts = set(random.sample(range(L), nr_muts))
            mut_seq = "".join([rec_seq[i] if i not in muts else random.choice([reverse_complement(rec_seq[i])]) for i in range(L)])
        elif args.exp_type == 'all':
            muts = set(random.sample(range(L), nr_muts))
            mut_seq = "".join([rec_seq[i] if i not in muts else random.choice(['', reverse_complement(rec_seq[i]), rec_seq[i] + random.choice("ACGT")]) for i in range(L)])
        else:
            print("Wrong --exp_type label specified")
            sys.exit()

        mut_rec_name = rec_acc + "|sim_copy_" + str(cnt) + "_nr_muts_" + str(nr_muts) + "|"
        cnt += 1
        new_rec[mut_rec_name] = mut_seq

    return new_rec




def main(args):
    mkdir_p(args.outfolder)
    orig_recs = { acc : seq for acc, (seq, _) in readfq(open(args.ref, 'r'))}

    # remove '-' in seqs if any
    for (acc, seq) in list(orig_recs.items()):
        orig_recs[acc] = "".join([s for s in seq if s != "-"])

    sim_recs = simulate_rec(args, orig_recs)

    outfile = open(os.path.join(args.outfolder, "sim_recombinants.fa"), "w")
    for acc, seq in sim_recs.items():
        outfile.write(">{0}\n{1}\n".format(acc, seq))

    sim_reads = simulate_reads(args, sim_recs)

    outfile = open(os.path.join(args.outfolder, "sim_reads.fq"), "w")
    for acc, (seq,qual) in sim_reads.items():
        outfile.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq, qual))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generates recombinants with SNPs at different abundancies and reads from them.")
    parser.add_argument('--ref', type=str, help='Path to fasta file with starting recombinans.')
    parser.add_argument('--C', type=int, default = 50, help='Total nr of output copies of recombinants (default 50).')
    parser.add_argument('--S', type=float, default = 0.01, help='Recombinant SNP/indel rate (default 0.01)')
    parser.add_argument('--N', type=int, default = 1000, help='Nr reads (default 1000)')
    parser.add_argument('--X', type=int, default = 4, help='Mean error rate (default 4). Choose error rate between 4, 7, or 12. ')
    parser.add_argument('--exp_type', type=str, default = 'subs', help='Experiment type. Choose "subs" for only SNP in Recombinants or "all" for SNP and indels (Default subs).')
    parser.add_argument('--outfolder', type=str, help='Outfolder.')
    
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)