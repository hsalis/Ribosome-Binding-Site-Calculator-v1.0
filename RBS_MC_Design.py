#This Python class generates a synthetic RBS sequence according to a user's target translation initiation rate, protein coding sequence, and pre-sequence. RBS constraints are not allowed.
#This file is part of the Ribosome Binding Site Calculator.

#The Ribosome Binding Site Calculator is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#The Ribosome Binding Site Calculator is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.

#Copyright 2008-2009 by Howard Salis.

from RBS_Calculator import RBS_Calculator
import random, math, sets

infinity = 1.0e20

#Constraints settings
max_kinetic_score = 0.50
max_three_state_indicator = 6.0
min_helical_loop = 4
max_helical_loop = 12
Max_RBS_Length = RBS_Calculator.cutoff

#dG ranges
dG_range_high = 25.0
dG_range_low = -18.0

num_rbs_calculations = 0

def dsu_sort(idx, seq):
    """Sorts a list of tuples according to the idx column using a Decorate-Sort-Undecorate method"""

    for i, e in enumerate(seq):
        seq[i] = (e[idx], e)
    seq.sort()
    seq.reverse()
    for i, e in enumerate(seq):
        seq[i] = e[1]
    return seq

def weighted_choice(list_of_tuples):
    """Randomly chooses from a list of choices according to their weighted probabilities."""

    n = random.uniform(0.0, 1.0)
    for item, weight in list_of_tuples:
        if n < weight:
            break
        n = n - weight
    return item

def Run_RBS_Calculator(pre_seq,post_seq,RBS,verbose=True):
    """Short cut function to run the RBS Calculator on a pre_sequence, CDS, and RBS."""

    if vars().has_key('estimator'): del(estimator)

    start_range = [len(pre_seq) + len(RBS) - 2, len(pre_seq) + len(RBS) + 2]

    mRNA = pre_seq.upper() + RBS + post_seq.upper()

    estimator = RBS_Calculator(mRNA, start_range, "")
    estimator.calc_dG()
    if verbose: estimator.print_dG()

    global num_rbs_calculations
    num_rbs_calculations+=1

    return estimator

def Generate_Random_RBS(All_Random = False, Max_length = 20, Pre_length = 5, PchooseSD = 0.5, Core_length = 6, max_nonoptimal_spacing = 5):
    """Generates a random RBS sequence tailored towards the target translation initiation rate."""

    RBS = []
    if All_Random:

        for i in range(Max_length):
            RBS.append(random.choice(["A","T","G","C"]))
        return "".join(RBS)

    for i in range(Pre_length):
        RBS.append(random.choice(["A","T","G","C"]))

    SD = ["T","A","A","G","G","A","G","G","T"]

    #Choose Core_length nucleotides. Choose from the SD sequence with probability PchooseSD
    #Choose from non-SD sequence with probability (1 - PchooseSD) / 3
    #The beginning/end of the Core_length wrt to the SD sequence is uniformly randomly determined

    Core_length = min(len(SD),Core_length) #Can't be greater then SD length
    diff = len(SD) - Core_length
    begin = int(random.random() * diff)

    for i in range(Core_length):
        rand = random.random()
        if rand <= PchooseSD:
            RBS.append(SD[begin+i])
        else:
            choices = ["A","T","G","C"]
            choices.remove(SD[begin+i])
            RBS.append(random.choice(choices))

    optimal_spacing = RBS_Calculator.optimal_spacing
    offset = diff - begin

    spacing = random.choice(range(max(0,offset + optimal_spacing - max_nonoptimal_spacing), offset + optimal_spacing + max_nonoptimal_spacing))

    for i in range(spacing):
        RBS.append(random.choice(["A","T","G","C"]))

    if len(RBS) > Max_length:
         RBS = RBS[len(RBS)-Max_length:len(RBS)+1]

    return "".join(RBS)

def calc_constraints(RBS,estimator):
    """Calculates the sequence constraints. Returns True if one is violated. """

    kinetic_score = estimator.kinetic_score_list[0]
    three_state_indicator = estimator.three_state_indicator_list[0]
    #(helical_loop_list,bulge_loop_list,helical_start_ends,bulge_start_ends) = estimator.calc_longest_loop_bulge(estimator,True,True,RBS)

    #print "KS = ", kinetic_score
    #print "3-state = ", three_state_indicator
    #print "max/min helical = ", max(helical_loop_list), min(helical_loop_list)

    if kinetic_score > max_kinetic_score: return True
    if three_state_indicator > max_three_state_indicator: return True
    #if min(helical_loop_list) < min_helical_loop: return True
    #if max(helical_loop_list) > max_helical_loop: return True

    return False

def RemoveStartCodons(sequence):
    """Removes any start codons from an input sequence."""

    import random
    import re

    regexp_str = "|".join(RBS_Calculator.start_codons)
    find_starts = re.compile(regexp_str)
    matches = find_starts.finditer(sequence.upper())

    new_seq = sequence[:]
    for match in matches:
        start_pos = match.start()

        triplet = []
        triplet.append(random.choice(['A','T','G','C']))
        triplet.append(random.choice(['A','G','C']))
        triplet.append(random.choice(['A','T','C']))

        new_seq = new_seq[0:start_pos] + "".join(triplet) + new_seq[start_pos+3:len(new_seq)+1]

    matches = find_starts.search(new_seq.upper())
    if matches is None:
        return new_seq
    else:
        return RemoveStartCodons(new_seq)

def compnt(nt):
    if (nt.upper() == 'A'): return 'T'
    if (nt.upper() == 'T'): return 'A'
    if (nt.upper() == 'G'): return 'C'
    if (nt.upper() == 'C'): return 'G'

def MCmove_lower_kinetic_score(pre_seq,post_seq,RBS,estimator = None,MaxIters=infinity):
    """Removes long-range base paired nucleotides from an mRNA sequence (pre-seq,CDS,RBS group). Used when generating initial conditions for synthetic RBS sequences."""

    if estimator is None:
        estimator = Run_RBS_Calculator(pre_seq,post_seq,RBS,verbose=False)

    counter=0
    kinetic_score = estimator.kinetic_score_list[0]
    while (counter < MaxIters and kinetic_score > max_kinetic_score):

        structure = estimator.mRNA_structure_list[0]
        mRNA = structure["mRNA"]
        RBS_begin = mRNA.find(RBS)
        RBS_end = RBS_begin + len(RBS)

        #Alter RBS to reduce kinetic score
        #Create a sorted list of bp'd nucleotides with the ones with the highest kinetic score at the top
        ks_list = []
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]
        for (nt_x,nt_y) in zip(bp_x,bp_y):
            ks_list.append( (nt_y - nt_x, nt_x, nt_y) )

        dsu_sort(0,ks_list) #Sort max-top according to 1st column: kinetic_score

        #Determine the bp'd nucleotides with the highest kinetic score. Are either located in the RBS?
        #If so, replace them with another random nucleotide
        nucleotides = sets.Set(['A', 'T', 'G', 'C'])

        num_mutations = min(len(ks_list),10)
        for i in range(num_mutations):
            nt_x = ks_list[i][1] - 1 #python index
            nt_y = ks_list[i][2] - 1 #python index

            #nt_x is in the RBS
            if (nt_x >= RBS_begin and nt_x < RBS_end):

                pos = nt_x - RBS_begin
                letter = random.choice(list(nucleotides ^ sets.Set(RBS[pos])))

                #print "Mutating ", RBS[pos], " --> ", letter

                RBS = RBS[0:pos] + letter + RBS[pos+1:len(RBS)+1]


            #nt_y is in the RBS
            elif (nt_y >= RBS_begin and nt_y < RBS_end):
                pos = nt_y - RBS_begin
                letter = random.choice(list(nucleotides ^ sets.Set(RBS[pos])))

                #print "Mutating ", RBS[pos], " --> ", letter

                RBS = RBS[0:pos] + letter + RBS[pos+1:len(RBS)+1]


            elif len(RBS) < RBS_Calculator.cutoff:
                #Insert a nucleotide at the 5' end of the RBS
                letter = random.choice(list(nucleotides))

                #print "Inserting ", letter, " at 5' end"

                RBS = letter + RBS

        RBS = RemoveStartCodons(RBS) #No start codons in RBS!!

        #print "RBS = ", RBS
        estimator = Run_RBS_Calculator(pre_seq,post_seq,RBS,verbose=False)
        kinetic_score = estimator.kinetic_score_list[0]

    return (RBS, estimator)


def MCmove_constrain_helical_loop(pre_seq,post_seq,RBS,estimator):
    """Modifies the mRNA sequence to reduce the size of helical loops."""

    #Alter RBS sequence so that max/min helical loop constraints are valid

    structure = estimator.mRNA_structure_list[0]
    (helical_loop_list,bulge_loop_list,helical_start_ends,bulge_start_ends) = estimator.calc_longest_loop_bulge(structure,True,True,RBS)

    RBS_begin = len(pre_seq)
    RBS_end = RBS_begin + len(RBS)

    #Insert or delete nucleotides to increase/decrease loop length
    for (loop_length,start_end) in zip(helical_loop_list,helical_start_ends):

        if loop_length > max_helical_loop:
            #Choose random nucleotide in loop. Delete it.

            #Identify what part of the loop is in the RBS
            RBS_range = sets.Set(range(RBS_begin+1,RBS_end+1))
            loop_range = sets.Set(range(start_end[0]+1,start_end[1]))
            change_range = list(RBS_range & loop_range) #Intersection

            #print "Loops in RBS: ", change_range

            if len(change_range) > 0:
                pos = random.choice(change_range) - len(pre_seq)
                RBS = RBS[0:pos] + RBS[pos+1:len(RBS)+1] #Delete nucleotide at position pos


        elif loop_length < min_helical_loop:
            #Choose random position in loop and insert a nucleotide before it.
            #Identify what part of the loop is in the RBS
            RBS_range = sets.Set(range(RBS_begin+1,RBS_end+1))
            loop_range = sets.Set(range(start_end[0]+1,start_end[1]))
            change_range = list(RBS_range & loop_range) #Intersection

            #print "Loops in RBS: ", change_range

            if len(change_range) > 0:
                pos = random.choice(change_range) - len(pre_seq)
                letter = random.choice(['A','T','C','G'])
                RBS = RBS[0:pos] + letter + RBS[pos+1:len(RBS)+1]

    #estimator = Run_RBS_Calculator(pre_seq,post_seq,RBS,verbose=False)

    return RBS

def GetInitialRBS(pre_seq,post_seq,dG_target):
    """Generates a random initial condition for designing a synthetic RBS sequence."""

    use_new = True
    Pre_length = 25

    dG_target_nondim = (dG_target - dG_range_high) / (dG_range_low - dG_range_high)
    #0.0: Low expression
    #1.0: High expression

    if     dG_target_nondim <  0.125:
        PchooseSD = 0.50
        Core_length = 4
        max_nonoptimal_spacing = 10
    elif   dG_target_nondim <  0.250:
        PchooseSD = 0.50
        Core_length = 4
        max_nonoptimal_spacing = 10
    elif   dG_target_nondim < 0.5:
        PchooseSD = 0.75
        Core_length = 4
        max_nonoptimal_spacing = 10
    elif   dG_target_nondim < 0.7:
        PchooseSD = 0.75
        Core_length = 4
        max_nonoptimal_spacing = 5
    elif   dG_target_nondim <  0.8:
        PchooseSD = 0.75
        Core_length = 6
        max_nonoptimal_spacing = 5
    elif   dG_target_nondim <  0.9:
        PchooseSD = 0.90
        Core_length = 6
        max_nonoptimal_spacing = 5
    elif   dG_target_nondim <  0.95:
        PchooseSD = 0.90
        Core_length = 8
        max_nonoptimal_spacing = 3
    elif   dG_target_nondim <=  1.00:
        PchooseSD = 1.0
        Core_length = 9
        max_nonoptimal_spacing = 2

    dG_total = dG_range_high + 1
    while dG_total > dG_range_high:
        RBS = Generate_Random_RBS(False, Max_RBS_Length, Pre_length, PchooseSD, Core_length, max_nonoptimal_spacing)
        RBS = RemoveStartCodons(RBS)

        estimator = Run_RBS_Calculator(pre_seq,post_seq,RBS,verbose=False)

        if use_new:
            RBS = MCmove_constrain_helical_loop(pre_seq,post_seq,RBS,estimator)
            (RBS,estimator) = MCmove_lower_kinetic_score(pre_seq,post_seq,RBS,estimator)

        else:
            counter=0
            while dG_total > 0 or calc_constraints(RBS,estimator):
                counter+=1
                RBS = Generate_Random_RBS(False, Max_RBS_Length, Pre_length, PchooseSD, Core_length, max_nonoptimal_spacing)
                RBS = RemoveStartCodons(RBS)

                estimator = Run_RBS_Calculator(pre_seq,post_seq,RBS,verbose=False)
                dG_total = estimator.dG_total_list[0]

        dG_total = estimator.dG_total_list[0]

    return (RBS,estimator)

def Monte_Carlo_Design(pre_seq, post_seq, RBS_init = None, TIR_target = None, dG_target = None, MaxIter = 10000, verbose = False):
    """Master function for designing synthetic RBS sequences without constraints."""

    #Check if dG_total or TIR (translation initiation rate) was specified. If TIR, then convert to dG_total.
    if TIR_target is not None:
        dG_target = RBS_Calculator.RT_eff * (RBS_Calculator.logK - math.log(float(TIR_target)))

    if verbose: print "dG_target = ", dG_target

    #Parameters
    max_init_energy = 10.0 #kcal/mol
    tol = 0.25 #kcal/mol
    annealing_accept_ratios = [0.01, 0.20] #first is min, second is max
    annealing_min_moves = 50
    RT_init = 0.6 #roughly 300K

    weighted_moves = [('insert',0.10),('delete',0.10),('replace',0.80)]

    #Define the energy/cost function based on the dG_target and the other, optional targets
    calc_energy = lambda (dG_total): abs(dG_total - dG_target)

    #If RBS_Init is given, use it. Otherwise, randomly choose one that is a decent starting point.
    if verbose: print "Determining Initial RBS"

    if RBS_init is None:
        (RBS,estimator) = GetInitialRBS(pre_seq,post_seq,dG_target)
    else:
        RBS = RBS_init
        estimator = Run_RBS_Calculator(pre_seq,post_seq,RBS,verbose=False)

    #Initialization
    counter = 0
    accepts = 0
    rejects = 0
    RT = RT_init

    dG_total = estimator.dG_total_list[0]
    energy = calc_energy(dG_total)

    if verbose: print "Initial RBS = ", RBS, " Energy = ", energy
    if verbose: estimator.print_dG(estimator.infinity)

    while energy > tol and counter < MaxIter:

        try:
            counter += 1
            accepted = False

            move = weighted_choice(weighted_moves)

            RBS_new = ''
            if verbose: print "Move #", counter, ": ", move
            if move == 'insert':

                pos = int(random.uniform(0.0,1.0) * len(RBS))
                letter = random.choice(['A', 'T', 'G', 'C'])
                RBS_new = RBS[0:pos] + letter + RBS[pos:len(RBS)]

            if move == 'delete':
                if (len(RBS) > 1):
                    pos = int(random.uniform(0.0,1.0) * len(RBS))
                    RBS_new = RBS[0:pos] + RBS[pos+1:len(RBS)]

            if move == 'replace':
                pos = int(random.uniform(0.0,1.0) * len(RBS))
                letter = random.choice(['A', 'T', 'G', 'C'])
                RBS_new = RBS[0:pos] + letter + RBS[pos+1:len(RBS)]

            RBS_new = RemoveStartCodons(RBS_new)

            if len(RBS_new) > Max_RBS_Length:
                RBS_new = RBS_new[len(RBS_new)-Max_RBS_Length:len(RBS_new)+1]

            estimator = Run_RBS_Calculator(pre_seq,post_seq,RBS_new,verbose=False)

            dG_total = estimator.dG_total_list[0]
            energy_new = calc_energy(dG_total)
            if calc_constraints(RBS_new,estimator):
                energy_new = infinity

            if verbose: print "New energy = ", energy_new

            if energy_new < energy:
                #Accept move immediately
                RBS = RBS_new
                energy = energy_new
                accepted = True
                if verbose: print "Move immediately accepted"
            else:

                ddE = (energy - energy_new)
                Metropolis = math.exp(ddE / RT)
                prob = random.uniform(0.0,1.0)

                if Metropolis > prob:
                    #Accept move based on conditional probability
                    RBS = RBS_new
                    energy = energy_new
                    accepts += 1
                    accepted = True
                    if verbose: print "Move conditionally accepted"
                else:
                    #Reject move
                    rejects += 1
                    if verbose: print "Move rejected"

            if accepted and verbose: estimator.print_dG(estimator.infinity)

            #Simulated annealing control
            if accepts + rejects > annealing_min_moves:

                ratio = float(accepts) / float(accepts + rejects)

                if ratio > annealing_accept_ratios[1]:
                    #Too many accepts, reduce RT
                    RT = RT / 2.0
                    accepts = 0
                    rejects = 0
                    if verbose: print "Accepting too many conditional moves, reducing temperature"

                if ratio < annealing_accept_ratios[0]:
                    #Too many rejects, increase RT
                    RT = RT * 2.0
                    accepts = 0
                    rejects = 0
                    if verbose: print "Rejecting too many conditional moves, increasing temperature"

        except KeyboardInterrupt:
            if verbose: print "Calculating Final State"
            estimator = Run_RBS_Calculator(pre_seq,post_seq,RBS,verbose=False)

            dG_total = estimator.dG_total_list[0]
            return (dG_total, RBS, estimator, counter)

    if verbose: estimator.print_dG(estimator.infinity)

    if verbose: print "Total number of RBS Evaluations: ", num_rbs_calculations

    if TIR_target is not None:

        #Convert back to TIR (Translation Initiation Rate)
        TIR_out = RBS_Calculator.K * math.exp(-dG_total / RBS_Calculator.RT_eff)
        return (TIR_out, RBS, estimator, counter)
    else:
        return (dG_total, RBS, estimator, counter)

def MC_Design_from_file(handle, output, dG_target, verbose = True, **kvars):
    """This function accepts a FASTA formatted file of Pre-sequences and CDSs and generates synthetic RBS sequences with the selected target dG_total. Uses Biopython for reading FASTA files. """

    from Bio import SeqIO

    #Set Defaults
    if not "kinetic_score_max" in kvars.keys(): kvars["kinetic_score_max"] = 0.500

    records = SeqIO.parse(handle,"fasta")

    for record in records:
        record_list.append(record)

    First = True
    counter=0

    assert len(dG_target) == len(record_list), "The length of dG_target should equal the number of designed RBSs in the input file"

    while counter < len(record_list):
        pre_seq = record_list[counter].seq.tostring().upper()
        post_seq = record_list[counter+1].seq.tostring().upper()
        counter += 2

        [dG_total, RBS, estimator] = Monte_Carlo_Design(pre_seq, post_seq, dG_target[counter/2-1], verbose,kvars)

        if verbose: estimator.print_dG(estimator.infinity)
        if verbose: print "Final RBS = ", RBS

        estimator.save_data(output, First)
        if First:
            First = False

        index = estimator.mRNA_rRNA_corrected_structure_list[0]["MinStructureID"]
        estimator.mRNA_rRNA_corrected_structure_list[0].export_PDF(index, "Monte Carlo Result", "MC_" + str(counter/2-1) + "_rRNA" + ".pdf")

        estimator.mRNA_structure_list[0].export_PDF(0, "Monte Carlo Result", "MC_" + str(counter/2-1) + "_mRNA" + ".pdf")

#-------------------------------------------------------------------------------
if __name__ == "__main__":


    pre_RBS = "TTCTAGA"
    post_RBS = "ATGCAGCACGTGTGCAGCACTACAGCGTGTGACGACTACAGCATTCACGACAGTCACATGCAGTTGACAC"
    dG_target = -10.0

    (dG_total, RBS, estimator, iterations) = Monte_Carlo_Design(pre_RBS, post_RBS, RBS_init = None, dG_target = dG_target, MaxIter = 10000, verbose = True)

    print "Finished"
    print "dG_total = ", dG_total
    print "RBS = ", RBS
    print "# iterations = ", iterations


