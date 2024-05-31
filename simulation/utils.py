def get_introgressed_tracts(ts, source_id, target_id):
    """
    Description:
        Outputs true introgressed tracts from a tree-sequence into a BED file.

    Arguments:
        ts tskit.TreeSequence: Tree-sequence containing introgressed tracts.
        chr_name int: Name of the chromosome.
        src_name str: Name of the source population.
        tgt_name str: Name of the target population.
        output string: Name of the output file.
    """
    introgressed_tracts = []
    for m in ts.migrations():
#         print(m)
        if m.dest == source_id and m.source == target_id: introgressed_tracts.append((int(m.left), int(m.right)))

    introgressed_tracts.sort(key=lambda x:x[0])
    return  introgressed_tracts

def getGraph_growth():
    """
    Description:
        Outputs demographic model sampled from unified distribution model with growth.
    """
    b = demes.Builder(
    description="Distribution of published demographic models"
    )
    
    if round(uniform(0,100))<75:
        CEU_anc = "OUT"
        ANC_end_time = round(uniform(2265, 15000))
        CEU_start_time = round(uniform(2265,ANC_end_time-1))
    else:
        CEU_anc = "ANC"
        CEU_start_time = round(uniform(2000,3000))
        ANC_end_time = CEU_start_time
    # ranges    
    pulse_time = round(uniform(1000,2000)) 
    pulse_proportions = round(uniform(0.02, 0.05), 2)
    ANC_start_size = round(uniform(10000, 25000))
    NEA_start_time = round(uniform(max(10000, ANC_end_time), 25000))
    NEA_start_size = round(uniform(2000, 10000))
    CEU_start_size = round(uniform(1000, 3500)) 
    CEU_exp_start_time = pulse_time - round(uniform(0,1000))
    CEU_exp_size = round(CEU_start_size*math.exp(CEU_exp_start_time*uniform(0.001, 0.005)))
    CEU_bn_start = round(uniform(pulse_time+40,  CEU_start_time))
    CEU_bn_end = CEU_bn_start - round(uniform(40,100))
    CEU_bn_size = round(uniform(100,1000))
    YRI_start_size = round(uniform(1000,10000))
    YRI_exp_start_time = round(uniform(0, 1000))
    YRI_exp_size =  round(YRI_start_size*math.exp(YRI_exp_start_time*uniform(0.001, 0.005)))
    b.add_deme("ANC", epochs=[dict(end_time=ANC_end_time, start_size=ANC_start_size)])
    b.add_deme("OUT", ancestors=["ANC"], epochs=[dict(start_size=YRI_start_size, end_time = YRI_exp_start_time),
                                                dict(start_size = YRI_start_size, end_size=YRI_exp_size,  end_time =0)])
    b.add_deme("SRC", ancestors=["ANC"], start_time=NEA_start_time, epochs=[dict(start_size=NEA_start_size)])
    b.add_deme("TAR", ancestors=[CEU_anc], start_time = CEU_start_time, epochs=[dict(end_time = CEU_bn_start, start_size=CEU_start_size), 
                                                                                dict(end_time = CEU_bn_end, start_size=CEU_bn_size),  
                                                                                dict(end_time = CEU_exp_start_time, start_size=CEU_start_size),
                                                                                dict(start_size = CEU_start_size, end_size=CEU_exp_size,  end_time =0),
                                                                               ])
    b.add_pulse(sources=["SRC"], dest="TAR", time=pulse_time, proportions= [ pulse_proportions])
    return b.resolve(), CEU_exp_size+YRI_exp_size

def getGraph_noGrowth():
    """
    Description:
        Outputs demographic model sampled from unified distribution model with no growth.
    """
    b = demes.Builder(
    description="Distribution of published demographic models"
    )
    
    if round(uniform(0,100))<75:
        CEU_anc = "OUT"
        ANC_end_time = round(uniform(2265, 15000))
        CEU_start_time = round(uniform(2265,ANC_end_time-1))
#         print("YRI")
    else:
        CEU_anc = "ANC"
        CEU_start_time = round(uniform(2000,3000))
        ANC_end_time = CEU_start_time
    # ranges    
    pulse_time = round(uniform(1000,2000)) 
    pulse_proportions = round(uniform(0.02, 0.05), 2)
    ANC_start_size = round(uniform(10000, 25000))
    NEA_start_time = round(uniform(max(10000, ANC_end_time), 25000))
    NEA_start_size = round(uniform(2000, 10000))
    CEU_start_size = round(uniform(1000, 10000)) 
    YRI_start_size = round(uniform(10000, 35000))
    b.add_deme("ANC", epochs=[dict(end_time=ANC_end_time, start_size=ANC_start_size)])
    b.add_deme("OUT", ancestors=["ANC"], epochs=[dict(start_size=YRI_start_size)])
    b.add_deme("SRC", ancestors=["ANC"], start_time=NEA_start_time, epochs=[dict(start_size=NEA_start_size)])
    b.add_deme("TAR", ancestors=[CEU_anc], start_time = CEU_start_time, epochs=[dict(start_size=CEU_start_size)])
    b.add_pulse(sources=["SRC"], dest="TAR", time=pulse_time, proportions= [ pulse_proportions])
    return b.resolve(), CEU_start_size+YRI_start_size